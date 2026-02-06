#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultN = 100'000'000'000ULL;
constexpr u64 kDefaultValidationN = 2'000ULL;
constexpr u64 kDirectHarmonicLimit = 5'000'000ULL;

// From the geometry of the problem (t = a / b), one gets:
// r(a,b) = b * (t - 2) / (t - 1) = b * (a - 2b) / (a - b).
long double r_formula(u64 a, u64 b) {
    return static_cast<long double>(b) * static_cast<long double>(a - 2 * b) /
           static_cast<long double>(a - b);
}

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
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

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2'000'000ULL) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    return std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

long double harmonic_asymptotic(u64 n) {
    if (n == 0) {
        return 0.0L;
    }

    // H_n = log(n) + gamma + 1/(2n) - 1/(12n^2) + 1/(120n^4) - ...
    constexpr long double kEulerGamma = 0.577215664901532860606512090082402431L;

    const long double x = static_cast<long double>(n);
    const long double inv = 1.0L / x;
    const long double inv2 = inv * inv;
    const long double inv4 = inv2 * inv2;
    const long double inv6 = inv4 * inv2;
    const long double inv8 = inv4 * inv4;
    const long double inv10 = inv8 * inv2;

    return std::log(x) + kEulerGamma + inv / 2.0L - inv2 / 12.0L + inv4 / 120.0L -
           inv6 / 252.0L + inv8 / 240.0L - 5.0L * inv10 / 660.0L;
}

long double harmonic_range(u64 left, u64 right) {
    if (left > right) {
        return 0.0L;
    }

    const u64 width = right - left + 1ULL;
    if (right <= kDirectHarmonicLimit || width <= kDirectHarmonicLimit) {
        long double sum = 0.0L;
        for (u64 k = left; k <= right; ++k) {
            sum += 1.0L / static_cast<long double>(k);
        }
        return sum;
    }

    return harmonic_asymptotic(right) - harmonic_asymptotic(left - 1ULL);
}

long double compute_g_fast(u64 n) {
    if (n < 3) {
        return 0.0L;
    }

    const u64 q = n / 2ULL;
    const u64 m = (n - 1ULL) / 2ULL;
    const u64 left = n - m;
    const u64 right = n - 1ULL;

    const long double qd = static_cast<long double>(q);
    const long double md = static_cast<long double>(m);
    const long double nd = static_cast<long double>(n);
    const long double ld = static_cast<long double>(left);
    const long double rd = static_cast<long double>(right);

    // Split by c = a - b:
    // c <= floor(n/2): closed polynomial contribution.
    const long double part_low =
        qd * (2.0L * qd * qd + 3.0L * qd - 5.0L) / 36.0L;

    // c > floor(n/2), with x = n - c in [1..m].
    const long double part_high_poly = md * (md + 1.0L) * (md + 2.0L) / 6.0L;

    const long double harmonic = harmonic_range(left, right);

    const long double acoef = nd * (2.0L * nd + 1.0L) * (nd + 1.0L);
    const long double bcoef = -(6.0L * nd * nd + 6.0L * nd + 1.0L);
    const long double ccoef = 6.0L * nd + 3.0L;

    const long double sum_k = (ld + rd) * md / 2.0L;
    const long double sum_k2 =
        (rd * (rd + 1.0L) * (2.0L * rd + 1.0L) -
         (ld - 1.0L) * ld * (2.0L * ld - 1.0L)) /
        6.0L;

    const long double rational_tail =
        acoef * harmonic + bcoef * md + ccoef * sum_k - 2.0L * sum_k2;

    return part_low + part_high_poly - rational_tail / 6.0L;
}

long double compute_g_direct_parallel(u64 n,
                                      bool allow_multithreading,
                                      unsigned requested_threads) {
    if (n < 3) {
        return 0.0L;
    }

    const long double pair_estimate =
        static_cast<long double>(n) * static_cast<long double>(n) / 4.0L;
    const std::size_t workload =
        (pair_estimate >= static_cast<long double>(std::numeric_limits<std::size_t>::max()))
            ? std::numeric_limits<std::size_t>::max()
            : static_cast<std::size_t>(pair_estimate);

    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, workload);
    if (threads == 1) {
        long double total = 0.0L;
        for (u64 a = 3; a <= n; ++a) {
            const u64 bmax = (a - 1ULL) / 2ULL;
            for (u64 b = 1; b <= bmax; ++b) {
                total += r_formula(a, b);
            }
        }
        return total;
    }

    const u64 a_count = n - 2ULL;
    std::vector<std::thread> workers;
    std::vector<long double> partial(threads, 0.0L);
    workers.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        const u64 a_begin = 3ULL + static_cast<u64>((static_cast<u128>(a_count) * t) / threads);
        const u64 a_end =
            3ULL + static_cast<u64>((static_cast<u128>(a_count) * (t + 1ULL)) / threads);

        workers.emplace_back([&, t, a_begin, a_end]() {
            long double local = 0.0L;
            for (u64 a = a_begin; a < a_end; ++a) {
                const u64 bmax = (a - 1ULL) / 2ULL;
                for (u64 b = 1; b <= bmax; ++b) {
                    local += r_formula(a, b);
                }
            }
            partial[t] = local;
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    long double total = 0.0L;
    for (const long double v : partial) {
        total += v;
    }
    return total;
}

bool nearly_equal(long double lhs,
                  long double rhs,
                  long double rel_tol = 1e-12L,
                  long double abs_tol = 1e-12L) {
    const long double diff = std::fabsl(lhs - rhs);
    if (diff <= abs_tol) {
        return true;
    }
    const long double scale = std::max(std::fabsl(lhs), std::fabsl(rhs));
    return diff <= rel_tol * scale;
}

std::string to_scientific_10_significant(long double value) {
    if (value == 0.0L) {
        return "0.000000000e0";
    }

    const bool negative = (value < 0.0L);
    long double x = negative ? -value : value;

    int exponent = static_cast<int>(std::floor(std::log10(x)));
    long double mantissa = x / std::powl(10.0L, static_cast<long double>(exponent));

    long double scaled = std::floor(mantissa * 1'000'000'000.0L + 0.5L);
    if (scaled >= 10'000'000'000.0L) {
        scaled /= 10.0L;
        ++exponent;
    }

    const u64 digits = static_cast<u64>(scaled);
    const u64 lead = digits / 1'000'000'000ULL;
    const u64 frac = digits % 1'000'000'000ULL;

    std::ostringstream oss;
    if (negative) {
        oss << '-';
    }
    oss << lead << '.' << std::setw(9) << std::setfill('0') << frac << 'e' << exponent;
    return oss.str();
}

struct Options {
    u64 n = kDefaultN;
    bool run_validation = true;
    u64 validation_n = kDefaultValidationN;
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
    bool answer_only = false;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-validation") {
            options.run_validation = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (arg == "--answer-only") {
            options.answer_only = true;
            continue;
        }

        u64 parsed_u64 = 0;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            options.n = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--validation-n=", parsed_u64)) {
            options.validation_n = parsed_u64;
            continue;
        }

        unsigned parsed_threads = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_threads)) {
            options.requested_threads = parsed_threads;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

bool run_checkpoints(const Options& options) {
    bool ok = true;

    auto check = [&](const bool condition, const std::string& label) {
        std::cout << "[checkpoint] " << label << ": " << (condition ? "PASS" : "FAIL") << '\n';
        if (!condition) {
            ok = false;
        }
    };

    check(nearly_equal(r_formula(3, 1), 0.5L, 1e-15L, 1e-15L), "r(3,1) = 1/2");
    check(nearly_equal(r_formula(6, 2), 1.0L, 1e-15L, 1e-15L), "r(6,2) = 1");
    check(nearly_equal(r_formula(12, 3), 2.0L, 1e-15L, 1e-15L), "r(12,3) = 2");

    const std::string g10 = to_scientific_10_significant(compute_g_fast(10));
    const std::string g100 = to_scientific_10_significant(compute_g_fast(100));

    check(g10 == "2.059722222e1", "G(10) scientific form");
    check(g100 == "1.922360980e4", "G(100) scientific form");

    if (options.validation_n >= 3) {
        const long double fast = compute_g_fast(options.validation_n);
        const long double direct = compute_g_direct_parallel(
            options.validation_n, options.allow_multithreading, options.requested_threads);
        const bool consistent = nearly_equal(fast, direct, 1e-12L, 1e-9L);

        std::ostringstream label;
        label << "fast vs direct at n=" << options.validation_n;
        check(consistent, label.str());

        if (!consistent) {
            std::cout << std::setprecision(20);
            std::cout << "  fast   = " << fast << '\n';
            std::cout << "  direct = " << direct << '\n';
        }
    }

    return ok;
}

} // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_validation && !options.answer_only) {
        if (!run_checkpoints(options)) {
            std::cerr << "Validation failed.\n";
            return 1;
        }
    }

    const long double answer = compute_g_fast(options.n);
    const std::string answer_scientific = to_scientific_10_significant(answer);

    if (!options.answer_only) {
        std::cout << "[result] G(" << options.n << ") = " << answer_scientific << '\n';
    } else {
        std::cout << answer_scientific << '\n';
    }

    return 0;
}
