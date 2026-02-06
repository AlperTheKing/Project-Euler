#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;

constexpr long double kDefaultA = 0.5L;
constexpr long double kDefaultTMax = 2000.0L;
constexpr long double kDefaultStep = 0.04L;
constexpr int kDefaultCutoff = 10'000;

constexpr long double kCoarseTMax = 1200.0L;
constexpr long double kCoarseStep = 0.05L;
constexpr int kCoarseCutoff = 4'000;

constexpr long double kFineTMax = 2000.0L;
constexpr long double kFineStep = 0.04L;
constexpr int kFineCutoff = 10'000;
constexpr long double kCheckpointExpectedMain = 0.565654540708545L;

struct Options {
    long double a = kDefaultA;
    long double t_max = kDefaultTMax;
    long double step = kDefaultStep;
    int cutoff = kDefaultCutoff;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct Kernel {
    explicit Kernel(const int cutoff)
        : cutoff_(cutoff), inv_half_squares_(static_cast<std::size_t>(cutoff) + 1ULL, 0.0L) {
        for (int n = 1; n <= cutoff_; ++n) {
            const long double nn = static_cast<long double>(n);
            inv_half_squares_[static_cast<std::size_t>(n)] = 0.5L / (nn * nn);
        }
    }

    long double psi(const long double t) const {
        long double p = 1.0L;
        for (int n = 1; n <= cutoff_; ++n) {
            p *= cosl(t * inv_half_squares_[static_cast<std::size_t>(n)]);
        }
        return p;
    }

  private:
    int cutoff_;
    std::vector<long double> inv_half_squares_;
};

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload) {
    if (!allow_multithreading || workload < 2ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::uint64_t parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_long_double_after_prefix(const std::string& arg,
                                    const char* prefix,
                                    long double& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    char* end = nullptr;
    errno = 0;
    const long double parsed = std::strtold(tail.c_str(), &end);
    if (end == tail.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(parsed)) {
        return false;
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        long double parsed_ld = 0.0L;
        if (parse_long_double_after_prefix(arg, "--a=", parsed_ld)) {
            options.a = parsed_ld;
            continue;
        }
        if (parse_long_double_after_prefix(arg, "--t-max=", parsed_ld)) {
            options.t_max = parsed_ld;
            continue;
        }
        if (parse_long_double_after_prefix(arg, "--step=", parsed_ld)) {
            options.step = parsed_ld;
            continue;
        }

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--cutoff=", parsed_u32)) {
            if (parsed_u32 > static_cast<u32>(std::numeric_limits<int>::max())) {
                std::cerr << "--cutoff is too large.\n";
                return false;
            }
            options.cutoff = static_cast<int>(parsed_u32);
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.step <= 0.0L) {
        std::cerr << "--step must be positive.\n";
        return false;
    }
    if (options.t_max <= 0.0L) {
        std::cerr << "--t-max must be positive.\n";
        return false;
    }
    if (options.cutoff < 1) {
        std::cerr << "--cutoff must be at least 1.\n";
        return false;
    }

    return true;
}

long double simpson_integral_parallel(const long double x,
                                      const long double t_max,
                                      const long double step,
                                      const int cutoff,
                                      const bool allow_multithreading,
                                      const unsigned requested_threads) {
    std::size_t segments = static_cast<std::size_t>(llround(t_max / step));
    if ((segments & 1ULL) != 0ULL) {
        ++segments;
    }

    const std::size_t interior_count = (segments > 1ULL) ? (segments - 1ULL) : 0ULL;
    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, interior_count);

    const Kernel kernel(cutoff);

    auto eval = [&](const std::size_t i) {
        if (i == 0ULL) {
            return x;
        }
        const long double t = static_cast<long double>(i) * step;
        return (sinl(t * x) / t) * kernel.psi(t);
    };

    std::vector<long double> partials(threads, 0.0L);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned tid = 0U; tid < threads; ++tid) {
        const std::size_t begin_offset = (interior_count * tid) / threads;
        const std::size_t end_offset = (interior_count * (tid + 1U)) / threads;

        const std::size_t begin_index = 1ULL + begin_offset;
        const std::size_t end_index = 1ULL + end_offset;

        pool.emplace_back([&, tid, begin_index, end_index]() {
            long double local = 0.0L;
            for (std::size_t i = begin_index; i < end_index; ++i) {
                const long double fi = eval(i);
                local += ((i & 1ULL) != 0ULL) ? (4.0L * fi) : (2.0L * fi);
            }
            partials[tid] = local;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    long double total = eval(0ULL) + eval(segments);
    for (const long double part : partials) {
        total += part;
    }

    return total * step / 3.0L;
}

long double probability_greater_than(const long double a,
                                     const long double t_max,
                                     const long double step,
                                     const int cutoff,
                                     const bool allow_multithreading,
                                     const unsigned requested_threads) {
    const long double pi = acosl(-1.0L);
    const long double zeta2 = (pi * pi) / 6.0L;
    const long double mu = zeta2 / 2.0L;

    if (a < 0.0L) {
        return 1.0L;
    }
    if (a >= zeta2) {
        return 0.0L;
    }

    const long double x = mu - a;
    const long double integral =
        simpson_integral_parallel(x, t_max, step, cutoff, allow_multithreading, requested_threads);

    return 0.5L + integral / pi;
}

bool close_to(const long double lhs, const long double rhs, const long double tolerance) {
    return fabsl(lhs - rhs) <= tolerance;
}

bool run_checkpoints(const Options& options) {
    const long double pi = acosl(-1.0L);
    const long double zeta2 = (pi * pi) / 6.0L;
    const long double mu = zeta2 / 2.0L;

    const long double p_mu = probability_greater_than(mu,
                                                      kFineTMax,
                                                      kFineStep,
                                                      kFineCutoff,
                                                      options.allow_multithreading,
                                                      options.requested_threads);
    if (!close_to(p_mu, 0.5L, 2.0e-12L)) {
        std::cerr << "Checkpoint failed: p(mu) should be 0.5, got " << std::setprecision(18)
                  << static_cast<double>(p_mu) << "\n";
        return false;
    }

    const long double p_zeta2_minus_1 = probability_greater_than(zeta2 - 1.0L,
                                                                  kFineTMax,
                                                                  kFineStep,
                                                                  kFineCutoff,
                                                                  options.allow_multithreading,
                                                                  options.requested_threads);
    if (!close_to(p_zeta2_minus_1, 0.5L, 2.0e-12L)) {
        std::cerr << "Checkpoint failed: p(zeta(2)-1) should be 0.5, got "
                  << std::setprecision(18) << static_cast<double>(p_zeta2_minus_1) << "\n";
        return false;
    }

    const long double p_one = probability_greater_than(1.0L,
                                                        kFineTMax,
                                                        kFineStep,
                                                        kFineCutoff,
                                                        options.allow_multithreading,
                                                        options.requested_threads);
    if (!close_to(p_one, 0.5L, 2.0e-12L)) {
        std::cerr << "Checkpoint failed: p(1) should be 0.5, got " << std::setprecision(18)
                  << static_cast<double>(p_one) << "\n";
        return false;
    }

    const long double symmetry_shift = 0.2L;
    const long double p_left = probability_greater_than(mu - symmetry_shift,
                                                         kCoarseTMax,
                                                         kCoarseStep,
                                                         kCoarseCutoff,
                                                         options.allow_multithreading,
                                                         options.requested_threads);
    const long double p_right = probability_greater_than(mu + symmetry_shift,
                                                          kCoarseTMax,
                                                          kCoarseStep,
                                                          kCoarseCutoff,
                                                          options.allow_multithreading,
                                                          options.requested_threads);
    if (!close_to(p_left + p_right, 1.0L, 2.0e-11L)) {
        std::cerr << "Checkpoint failed: symmetry mismatch, p(mu-u)+p(mu+u)="
                  << std::setprecision(18) << static_cast<double>(p_left + p_right) << "\n";
        return false;
    }

    const long double coarse_main = probability_greater_than(0.5L,
                                                              kCoarseTMax,
                                                              kCoarseStep,
                                                              kCoarseCutoff,
                                                              options.allow_multithreading,
                                                              options.requested_threads);
    const long double fine_main = probability_greater_than(0.5L,
                                                            kFineTMax,
                                                            kFineStep,
                                                            kFineCutoff,
                                                            options.allow_multithreading,
                                                            options.requested_threads);

    if (!close_to(fine_main, kCheckpointExpectedMain, 5.0e-11L)) {
        std::cerr << "Checkpoint failed: p(0.5) reference mismatch, got "
                  << std::setprecision(18) << static_cast<double>(fine_main) << "\n";
        return false;
    }

    if (!close_to(coarse_main, fine_main, 2.0e-10L)) {
        std::cerr << "Checkpoint failed: coarse/fine mismatch, coarse=" << std::setprecision(18)
                  << static_cast<double>(coarse_main) << ", fine=" << static_cast<double>(fine_main)
                  << "\n";
        return false;
    }

    const unsigned thread_probe =
        choose_thread_count(true, options.requested_threads, static_cast<std::size_t>(50'000ULL));
    if (options.allow_multithreading && thread_probe > 1U) {
        const long double single_thread =
            probability_greater_than(0.5L, kCoarseTMax, kCoarseStep, kCoarseCutoff, false, 1U);
        const long double multi_thread = probability_greater_than(0.5L,
                                                                   kCoarseTMax,
                                                                   kCoarseStep,
                                                                   kCoarseCutoff,
                                                                   true,
                                                                   thread_probe);
        if (!close_to(single_thread, multi_thread, 1.0e-13L)) {
            std::cerr << "Checkpoint failed: single/multi thread mismatch, single="
                      << std::setprecision(18) << static_cast<double>(single_thread)
                      << ", multi=" << static_cast<double>(multi_thread) << "\n";
            return false;
        }
    }

    return true;
}

} // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const auto start = std::chrono::steady_clock::now();

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const long double probability = probability_greater_than(options.a,
                                                              options.t_max,
                                                              options.step,
                                                              options.cutoff,
                                                              options.allow_multithreading,
                                                              options.requested_threads);

    const auto finish = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(finish - start).count();

    if (options.run_checkpoints) {
        std::cout << "Checkpoints passed.\n";
    }

    std::cout << std::fixed << std::setprecision(12)
              << "p(" << static_cast<double>(options.a) << ") = "
              << static_cast<double>(probability) << "\n";

    std::cout << std::fixed << std::setprecision(8)
              << "Rounded to 8 decimals: " << static_cast<double>(probability) << "\n";

    std::cout << std::fixed << std::setprecision(6)
              << "Elapsed: " << elapsed << " s\n";

    return 0;
}
