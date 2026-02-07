#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kDefaultN = 10'000'000U;
constexpr u32 kParallelThreshold = 200'000U;
constexpr u32 kChunkSize = 2'048U;

struct Options {
    u32 n = kDefaultN;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
    bool run_checks = true;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value_out) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value_out = parsed;
    return true;
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value_out) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<u32>::max())) {
        return false;
    }
    value_out = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value_out) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value_out = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (arg == "--skip-checks") {
            options.run_checks = false;
            continue;
        }

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--n=", parsed_u32)) {
            options.n = parsed_u32;
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

    if (options.n < 2U) {
        std::cerr << "--n must be at least 2.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const u64 workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) {
        return 1U;
    }

    const unsigned max_threads = static_cast<unsigned>(
        std::min<u64>(workload_units, static_cast<u64>(std::numeric_limits<unsigned>::max())));

    if (requested_threads > 0U) {
        return std::max(1U, std::min(requested_threads, max_threads));
    }

    if (workload_units < static_cast<u64>(kParallelThreshold)) {
        return 1U;
    }

    unsigned hardware_threads = std::thread::hardware_concurrency();
    if (hardware_threads == 0U) {
        hardware_threads = 1U;
    }

    return std::max(1U, std::min(hardware_threads, max_threads));
}

struct KahanAccumulator {
    long double sum = 0.0L;
    long double compensation = 0.0L;

    void add(const long double value) {
        const long double y = value - compensation;
        const long double t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
    }
};

class Euler441Solver {
public:
    explicit Euler441Solver(const u32 n)
        : n_(n),
          half_(n / 2U),
          spf_(static_cast<std::size_t>(n) + 1U, 0U),
          phi_(static_cast<std::size_t>(n) + 1U, 0U),
          harmonic_(static_cast<std::size_t>(n) + 1U, 0.0L) {
        build_linear_sieve();
        build_harmonic_prefix();
    }

    long double solve(const unsigned thread_count) const {
        if (n_ < 2U) {
            return 0.0L;
        }

        const unsigned clamped_threads =
            std::max(1U, std::min(thread_count, static_cast<unsigned>(n_ - 1U)));

        if (clamped_threads == 1U) {
            KahanAccumulator acc;
            for (u32 q = 2U; q <= n_; ++q) {
                acc.add(compute_contribution(q));
            }
            return acc.sum;
        }

        std::atomic<u32> next_q(2U);
        std::vector<long double> partial(clamped_threads, 0.0L);
        std::vector<std::thread> workers;
        workers.reserve(clamped_threads);

        for (unsigned tid = 0U; tid < clamped_threads; ++tid) {
            workers.emplace_back([&, tid]() {
                KahanAccumulator local;

                while (true) {
                    const u32 q_begin = next_q.fetch_add(kChunkSize, std::memory_order_relaxed);
                    if (q_begin > n_) {
                        break;
                    }

                    const u32 q_end = std::min<u32>(n_, q_begin + kChunkSize - 1U);
                    for (u32 q = q_begin; q <= q_end; ++q) {
                        local.add(compute_contribution(q));
                    }
                }

                partial[static_cast<std::size_t>(tid)] = local.sum;
            });
        }

        for (auto& worker : workers) {
            worker.join();
        }

        KahanAccumulator total;
        for (const long double value : partial) {
            total.add(value);
        }
        return total.sum;
    }

private:
    static constexpr int kMaxPrimeFactors = 16;
    static constexpr int kMaxSquarefreeDivisors = 1024;

    u32 n_ = 0U;
    u32 half_ = 0U;
    std::vector<u32> spf_;
    std::vector<u32> phi_;
    std::vector<long double> harmonic_;

    void build_linear_sieve() {
        std::vector<u32> primes;
        primes.reserve(static_cast<std::size_t>(n_) / 10U + 16U);

        phi_[1] = 1U;

        for (u32 i = 2U; i <= n_; ++i) {
            if (spf_[i] == 0U) {
                spf_[i] = i;
                phi_[i] = i - 1U;
                primes.push_back(i);
            }

            for (const u32 p : primes) {
                const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
                if (v > n_) {
                    break;
                }

                spf_[static_cast<std::size_t>(v)] = p;
                if (p == spf_[i]) {
                    phi_[static_cast<std::size_t>(v)] = phi_[i] * p;
                    break;
                }

                phi_[static_cast<std::size_t>(v)] = phi_[i] * (p - 1U);
            }
        }
    }

    void build_harmonic_prefix() {
        harmonic_[0] = 0.0L;
        for (u32 i = 1U; i <= n_; ++i) {
            harmonic_[static_cast<std::size_t>(i)] =
                harmonic_[static_cast<std::size_t>(i - 1U)] + 1.0L / static_cast<long double>(i);
        }
    }

    long double compute_contribution(const u32 q) const {
        u32 prime_factors[kMaxPrimeFactors];
        int factor_count = 0;

        u32 x = q;
        while (x > 1U) {
            const u32 p = spf_[x];
            prime_factors[factor_count++] = p;
            while (x % p == 0U) {
                x /= p;
            }
        }

        u32 divisors[kMaxSquarefreeDivisors];
        int signs[kMaxSquarefreeDivisors];
        int subset_count = 1;

        divisors[0] = 1U;
        signs[0] = 1;

        for (int i = 0; i < factor_count; ++i) {
            const u32 p = prime_factors[i];
            for (int j = 0; j < subset_count; ++j) {
                divisors[subset_count + j] = divisors[j] * p;
                signs[subset_count + j] = -signs[j];
            }
            subset_count <<= 1;
        }

        const u32 qm1 = q - 1U;
        long double a_full = 0.0L;
        for (int i = 0; i < subset_count; ++i) {
            const u32 d = divisors[i];
            const long double signed_inv_d =
                static_cast<long double>(signs[i]) / static_cast<long double>(d);
            a_full += signed_inv_d * harmonic_[static_cast<std::size_t>(qm1 / d)];
        }

        if (q <= half_) {
            return (static_cast<long double>(phi_[q]) + a_full) / static_cast<long double>(q);
        }

        const u32 a = n_ - q;
        long double a_a = 0.0L;
        long double b = 0.0L;
        i64 c = 0;

        for (int i = 0; i < subset_count; ++i) {
            const u32 d = divisors[i];
            const int sign = signs[i];
            const long double signed_inv_d =
                static_cast<long double>(sign) / static_cast<long double>(d);
            const u32 a_over_d = a / d;

            c += static_cast<i64>(sign) * static_cast<i64>(a_over_d);
            const long double ha = harmonic_[static_cast<std::size_t>(a_over_d)];
            const long double hq = harmonic_[static_cast<std::size_t>(qm1 / d)];

            a_a += signed_inv_d * ha;
            b += signed_inv_d * (hq - ha);
        }

        const long double numerator =
            static_cast<long double>(c) + a_a + static_cast<long double>(a + 1U) * b;
        return numerator / static_cast<long double>(q);
    }
};

long double brute_force_s(const u32 n) {
    long double total = 0.0L;

    for (u32 i = 2U; i <= n; ++i) {
        long double r_i = 0.0L;
        for (u32 p = 1U; p < i; ++p) {
            for (u32 q = p + 1U; q <= i; ++q) {
                if (p + q < i) {
                    continue;
                }
                if (std::gcd(p, q) != 1U) {
                    continue;
                }
                r_i += 1.0L / (static_cast<long double>(p) * static_cast<long double>(q));
            }
        }
        total += r_i;
    }

    return total;
}

bool near_abs(const long double a, const long double b, const long double eps) {
    return std::fabsl(a - b) <= eps;
}

bool run_problem_checkpoints() {
    struct Checkpoint {
        u32 n;
        long double expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {2U, 0.5L},
        {10U, 6.9146825396825396825L},
        {100U, 58.296238062166323935L},
    };

    for (const Checkpoint& checkpoint : checkpoints) {
        const Euler441Solver solver(checkpoint.n);
        const long double got = solver.solve(1U);
        if (!near_abs(got, checkpoint.expected, 1e-12L)) {
            std::cerr << "Checkpoint failed for n=" << checkpoint.n << ": got "
                      << std::setprecision(18) << static_cast<double>(got) << ", expected "
                      << static_cast<double>(checkpoint.expected) << '\n';
            return false;
        }
    }

    return true;
}

bool run_bruteforce_crosscheck() {
    for (u32 n = 2U; n <= 40U; ++n) {
        const Euler441Solver solver(n);
        const long double fast = solver.solve(1U);
        const long double brute = brute_force_s(n);
        if (!near_abs(fast, brute, 1e-12L)) {
            std::cerr << "Bruteforce mismatch at n=" << n << ": fast=" << std::setprecision(18)
                      << static_cast<double>(fast) << ", brute=" << static_cast<double>(brute)
                      << '\n';
            return false;
        }
    }
    return true;
}

bool run_parallel_consistency_check() {
    unsigned hardware_threads = std::thread::hardware_concurrency();
    if (hardware_threads < 2U) {
        return true;
    }

    const u32 n = 200'000U;
    const Euler441Solver solver(n);
    const long double single = solver.solve(1U);

    const unsigned threads = std::max(2U, std::min(8U, hardware_threads));
    const long double parallel = solver.solve(threads);

    if (!near_abs(single, parallel, 1e-12L)) {
        std::cerr << "Parallel consistency check failed: single=" << std::setprecision(18)
                  << static_cast<double>(single) << ", parallel="
                  << static_cast<double>(parallel) << ", threads=" << threads << '\n';
        return false;
    }

    return true;
}

bool run_validation_suite(const bool allow_multithreading) {
    if (!run_problem_checkpoints()) {
        return false;
    }
    if (!run_bruteforce_crosscheck()) {
        return false;
    }
    if (allow_multithreading && !run_parallel_consistency_check()) {
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checks && !run_validation_suite(options.allow_multithreading)) {
        return 1;
    }

    const Euler441Solver solver(options.n);
    const unsigned threads = choose_thread_count(options.allow_multithreading,
                                                 options.requested_threads,
                                                 static_cast<u64>(options.n) - 1ULL);
    const long double answer = solver.solve(threads);

    std::cout << std::fixed << std::setprecision(4) << answer << '\n';
    return 0;
}
