#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 100'000'000ULL;
constexpr u32 kDefaultNMin = 1U;
constexpr u32 kDefaultNMax = 30U;
constexpr u64 kDefaultExponent = 987'654'321ULL;
constexpr std::size_t kAutoThreadThreshold = 96ULL;

struct Options {
    u32 n_min = kDefaultNMin;
    u32 n_max = kDefaultNMax;
    u64 exponent = kDefaultExponent;
    bool run_checks = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
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

    value = parsed;
    return true;
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<u32>::max())) {
        return false;
    }
    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg, const char* prefix, unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checks") {
            options.run_checks = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--n-min=", parsed_u32)) {
            options.n_min = parsed_u32;
            continue;
        }
        if (parse_u32_after_prefix(arg, "--n-max=", parsed_u32)) {
            options.n_max = parsed_u32;
            continue;
        }

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--exp=", parsed_u64)) {
            options.exponent = parsed_u64;
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

    if (options.n_min == 0U || options.n_max == 0U) {
        std::cerr << "n bounds must be >= 1.\n";
        return false;
    }
    if (options.n_min > options.n_max) {
        std::cerr << "Invalid range: --n-min must be <= --n-max.\n";
        return false;
    }
    if (options.exponent == 0ULL || (options.exponent & 1ULL) == 0ULL) {
        std::cerr << "Exponent must be a positive odd integer.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) {
        return 1U;
    }

    const unsigned max_threads = static_cast<unsigned>(
        std::min<std::size_t>(workload_units, static_cast<std::size_t>(std::numeric_limits<unsigned>::max())));

    if (requested_threads > 0U) {
        return std::max(1U, std::min(requested_threads, max_threads));
    }

    if (workload_units < kAutoThreadThreshold) {
        return 1U;
    }

    unsigned hardware_threads = std::thread::hardware_concurrency();
    if (hardware_threads == 0U) {
        hardware_threads = 1U;
    }
    return std::max(1U, std::min(hardware_threads, max_threads));
}

u64 pow_mod_u64(u64 base, u64 exponent, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }

    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    u64 exp = exponent;

    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        exp >>= 1ULL;
    }

    return result;
}

struct Matrix3 {
    u64 v[3][3] = {{0ULL, 0ULL, 0ULL}, {0ULL, 0ULL, 0ULL}, {0ULL, 0ULL, 0ULL}};
};

Matrix3 make_identity_matrix3() {
    Matrix3 out;
    for (int i = 0; i < 3; ++i) {
        out.v[i][i] = 1ULL;
    }
    return out;
}

Matrix3 multiply_matrix3(const Matrix3& a, const Matrix3& b, const u64 mod) {
    Matrix3 out;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            u128 sum = 0;
            for (int k = 0; k < 3; ++k) {
                sum += static_cast<u128>(a.v[i][k]) * static_cast<u128>(b.v[k][j]);
            }
            out.v[i][j] = static_cast<u64>(sum % mod);
        }
    }
    return out;
}

Matrix3 power_matrix3(Matrix3 base, u64 exponent, const u64 mod) {
    Matrix3 result = make_identity_matrix3();
    u64 exp = exponent;

    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = multiply_matrix3(result, base, mod);
        }
        base = multiply_matrix3(base, base, mod);
        exp >>= 1ULL;
    }

    return result;
}

u64 compute_newton_sum_mod(const u32 n, const u64 exponent, const u64 mod) {
    const u64 c = pow_mod_u64(2ULL, static_cast<u64>(n), mod);
    const u64 s0 = 3ULL % mod;
    const u64 s1 = c;
    const u64 s2 = static_cast<u64>((static_cast<u128>(c) * c) % mod);

    if (exponent == 0ULL) {
        return s0;
    }
    if (exponent == 1ULL) {
        return s1;
    }
    if (exponent == 2ULL) {
        return s2;
    }

    Matrix3 transition;
    transition.v[0][0] = c;
    transition.v[0][1] = 0ULL;
    transition.v[0][2] = (mod - (static_cast<u64>(n) % mod)) % mod;
    transition.v[1][0] = 1ULL % mod;
    transition.v[1][1] = 0ULL;
    transition.v[1][2] = 0ULL;
    transition.v[2][0] = 0ULL;
    transition.v[2][1] = 1ULL % mod;
    transition.v[2][2] = 0ULL;

    const Matrix3 p = power_matrix3(transition, exponent - 2ULL, mod);
    const u128 top = static_cast<u128>(p.v[0][0]) * s2 +
                     static_cast<u128>(p.v[0][1]) * s1 +
                     static_cast<u128>(p.v[0][2]) * s0;
    return static_cast<u64>(top % mod);
}

u64 compute_newton_sum_mod_iterative(const u32 n, const u64 exponent, const u64 mod) {
    const u64 c = pow_mod_u64(2ULL, static_cast<u64>(n), mod);
    const u64 n_mod = static_cast<u64>(n) % mod;
    const u64 s0 = 3ULL % mod;
    const u64 s1 = c;
    const u64 s2 = static_cast<u64>((static_cast<u128>(c) * c) % mod);

    if (exponent == 0ULL) {
        return s0;
    }
    if (exponent == 1ULL) {
        return s1;
    }
    if (exponent == 2ULL) {
        return s2;
    }

    u64 a = s0;
    u64 b = s1;
    u64 cur = s2;
    for (u64 k = 3ULL; k <= exponent; ++k) {
        const u64 lhs = static_cast<u64>((static_cast<u128>(c) * cur) % mod);
        const u64 rhs = static_cast<u64>((static_cast<u128>(n_mod) * a) % mod);
        u64 next = lhs;
        if (next >= rhs) {
            next -= rhs;
        } else {
            next += mod - rhs;
        }
        a = b;
        b = cur;
        cur = next;
    }
    return cur;
}

u128 compute_newton_sum_exact_small(const u32 n, const u32 exponent) {
    u128 c = 1;
    for (u32 i = 0; i < n; ++i) {
        c *= 2U;
    }

    const u128 s0 = 3U;
    const u128 s1 = c;
    const u128 s2 = c * c;

    if (exponent == 0U) {
        return s0;
    }
    if (exponent == 1U) {
        return s1;
    }
    if (exponent == 2U) {
        return s2;
    }

    u128 a = s0;
    u128 b = s1;
    u128 cur = s2;
    for (u32 k = 3U; k <= exponent; ++k) {
        const u128 next = c * cur - static_cast<u128>(n) * a;
        a = b;
        b = cur;
        cur = next;
    }
    return cur;
}

long double evaluate_polynomial(const long double x, const long double c, const u32 n) {
    return x * x * x - c * x * x + static_cast<long double>(n);
}

long double largest_real_root(const u32 n) {
    const long double c = std::ldexp(1.0L, static_cast<int>(n));
    long double low = std::max(0.0L, c - 1.0L);
    long double high = c;

    while (evaluate_polynomial(low, c, n) > 0.0L && low > 0.0L) {
        low = std::max(0.0L, low - 1.0L);
    }

    for (int iter = 0; iter < 220; ++iter) {
        const long double mid = (low + high) * 0.5L;
        if (evaluate_polynomial(mid, c, n) <= 0.0L) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return (low + high) * 0.5L;
}

u64 floor_power(const long double base, const u32 exponent) {
    long double value = 1.0L;
    for (u32 i = 0; i < exponent; ++i) {
        value *= base;
    }
    const long double floored = std::floor(value + 1e-18L);
    return (floored <= 0.0L) ? 0ULL : static_cast<u64>(floored);
}

u64 compute_term_mod(const u32 n, const u64 exponent) {
    const u64 sum_mod = compute_newton_sum_mod(n, exponent, kMod);
    return (sum_mod + kMod - 1ULL) % kMod;
}

u64 solve_range(const u32 n_min,
                const u32 n_max,
                const u64 exponent,
                const bool allow_multithreading,
                const unsigned requested_threads) {
    const std::size_t workload = static_cast<std::size_t>(n_max) - static_cast<std::size_t>(n_min) + 1ULL;
    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, workload);

    if (threads <= 1U) {
        u64 total = 0ULL;
        for (u32 n = n_min; n <= n_max; ++n) {
            total += compute_term_mod(n, exponent);
            if (total >= kMod) {
                total -= kMod;
            }
        }
        return total;
    }

    std::atomic<u32> next_n(n_min);
    std::vector<u64> partial(static_cast<std::size_t>(threads), 0ULL);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    auto worker = [&](const unsigned tid) {
        u64 local = 0ULL;
        while (true) {
            const u32 n = next_n.fetch_add(1U, std::memory_order_relaxed);
            if (n > n_max) {
                break;
            }
            local += compute_term_mod(n, exponent);
            if (local >= kMod) {
                local -= kMod;
            }
        }
        partial[static_cast<std::size_t>(tid)] = local;
    };

    for (unsigned t = 0U; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    for (std::thread& th : pool) {
        th.join();
    }

    u64 total = 0ULL;
    for (const u64 part : partial) {
        total += part;
        if (total >= kMod) {
            total -= kMod;
        }
    }
    return total;
}

u64 solve(const Options& options) {
    return solve_range(options.n_min,
                       options.n_max,
                       options.exponent,
                       options.allow_multithreading,
                       options.requested_threads);
}

bool run_validations(const Options& options) {
    // Checkpoint 1: statement example floor(a_2^2) = 14.
    {
        const long double a2 = largest_real_root(2U);
        const u64 floor_a2_sq = floor_power(a2, 2U);
        if (floor_a2_sq != 14ULL) {
            std::cerr << "Validation failed (1/4): floor(a_2^2) expected 14, got "
                      << floor_a2_sq << '\n';
            return false;
        }
        std::cout << "Validation 1/4 [PASS]\n";
    }

    // Checkpoint 2: matrix recurrence matches direct iterative recurrence.
    {
        for (u32 n = 1U; n <= 12U; ++n) {
            for (u64 k = 0ULL; k <= 220ULL; ++k) {
                const u64 by_matrix = compute_newton_sum_mod(n, k, kMod);
                const u64 by_iter = compute_newton_sum_mod_iterative(n, k, kMod);
                if (by_matrix != by_iter) {
                    std::cerr << "Validation failed (2/4): recurrence mismatch at n=" << n
                              << ", k=" << k << " (matrix=" << by_matrix
                              << ", iterative=" << by_iter << ")\n";
                    return false;
                }
            }
        }
        std::cout << "Validation 2/4 [PASS]\n";
    }

    // Checkpoint 3: floor(a_n^m) == S_m(n)-1 for small odd-exponent spot checks.
    {
        const u32 odd_exponents[] = {1U, 3U, 5U, 7U, 9U, 11U};
        for (u32 n = 1U; n <= 3U; ++n) {
            const long double alpha = largest_real_root(n);
            for (const u32 m : odd_exponents) {
                const u64 floor_numeric = floor_power(alpha, m);
                const u128 exact_sum = compute_newton_sum_exact_small(n, m);
                const u64 floor_via_identity = static_cast<u64>(exact_sum - 1U);
                if (floor_numeric != floor_via_identity) {
                    std::cerr << "Validation failed (3/4): floor identity mismatch at n=" << n
                              << ", m=" << m << " (numeric=" << floor_numeric
                              << ", identity=" << floor_via_identity << ")\n";
                    return false;
                }
            }
        }
        std::cout << "Validation 3/4 [PASS]\n";
    }

    // Checkpoint 4: thread consistency on a larger reduced workload.
    if (options.allow_multithreading) {
        unsigned hardware_threads = std::thread::hardware_concurrency();
        if (hardware_threads == 0U) {
            hardware_threads = 1U;
        }

        if (hardware_threads > 1U) {
            const u32 check_n_min = 1U;
            const u32 check_n_max = 400U;
            const u64 single = solve_range(check_n_min, check_n_max, options.exponent, false, 1U);
            const unsigned compare_threads =
                (options.requested_threads > 1U) ? options.requested_threads : 2U;
            const u64 multi =
                solve_range(check_n_min, check_n_max, options.exponent, true, compare_threads);
            if (single != multi) {
                std::cerr << "Validation failed (4/4): single-thread vs multi-thread mismatch"
                          << " (single=" << single << ", multi=" << multi << ")\n";
                return false;
            }
            std::cout << "Validation 4/4 [PASS]\n";
        } else {
            std::cout << "Validation 4/4 [SKIP: single-core environment]\n";
        }
    } else {
        std::cout << "Validation 4/4 [SKIP: --single-thread enabled]\n";
    }

    return true;
}

void print_usage() {
    std::cerr << "Usage: ./Euler356 [--skip-checks] [--single-thread] [--threads=<unsigned>] "
                 "[--n-min=<u32>] [--n-max=<u32>] [--exp=<odd-u64>]\n";
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        print_usage();
        return 1;
    }

    if (options.run_checks && !run_validations(options)) {
        return 1;
    }

    const u64 answer = solve(options);
    std::cout << std::setfill('0') << std::setw(8) << answer << '\n';
    return 0;
}
