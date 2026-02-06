#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr u64 kDefaultN = 1'000'000'000ULL;
constexpr int kModulus = 77'777'777;
constexpr int kPrimeFactors[5] = {7, 11, 73, 101, 137};

struct Options {
    u64 n = kDefaultN;
    bool run_checkpoints = true;
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

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
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

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            options.n = parsed_u64;
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

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

int normalize_mod(i64 x, const int mod) {
    const i64 r = x % static_cast<i64>(mod);
    return static_cast<int>(r < 0 ? r + mod : r);
}

int pow_mod(int base, int exp, const int mod) {
    i64 result = 1;
    i64 b = normalize_mod(base, mod);

    while (exp > 0) {
        if ((exp & 1) != 0) {
            result = (result * b) % mod;
        }
        b = (b * b) % mod;
        exp >>= 1;
    }

    return static_cast<int>(result);
}

int inverse_mod_prime(const int a, const int p) {
    // Fermat inverse, valid because p is prime and gcd(a, p)=1.
    return pow_mod(a, p - 2, p);
}

int factorial_mod_small(const u64 n, const int p) {
    int f = 1 % p;
    for (u64 i = 1ULL; i <= n; ++i) {
        f = static_cast<int>((static_cast<i64>(f) * static_cast<i64>(i % static_cast<u64>(p))) % p);
    }
    return f;
}

std::vector<int> compute_A_mod_prime(const int p, const u64 max_n) {
    std::vector<int> values(static_cast<std::size_t>(max_n) + 1ULL, 0);
    std::vector<int> comb(static_cast<std::size_t>(max_n) + 1ULL, 0);
    comb[0] = 1;

    int fact = 1 % p;

    for (u64 n = 0ULL; n <= max_n; ++n) {
        if (n > 0ULL) {
            comb[static_cast<std::size_t>(n)] = 1;
            for (u64 k = n - 1ULL; k >= 1ULL; --k) {
                const std::size_t idx = static_cast<std::size_t>(k);
                comb[idx] += comb[idx - 1ULL];
                if (comb[idx] >= p) {
                    comb[idx] -= p;
                }
            }

            if (n < static_cast<u64>(p)) {
                fact = static_cast<int>((static_cast<i64>(fact) * static_cast<i64>(n)) % p);
            } else {
                fact = 0;
            }
        }

        i64 current = fact;
        for (u64 k = 0ULL; k < n; ++k) {
            current += static_cast<i64>(comb[static_cast<std::size_t>(k)]) *
                       static_cast<i64>(values[static_cast<std::size_t>(k)]);
            current %= p;
        }

        values[static_cast<std::size_t>(n)] = static_cast<int>(current);
    }

    return values;
}

u64 reduced_index_for_prime(const u64 n, const int p) {
    if (n < static_cast<u64>(p)) {
        return n;
    }

    const u64 period = static_cast<u64>(p) * static_cast<u64>(p - 1);
    return static_cast<u64>(p) + ((n - static_cast<u64>(p)) % period);
}

int solve_mod_prime(const int p, const u64 n) {
    const u64 index = reduced_index_for_prime(n, p);
    const std::vector<int> A = compute_A_mod_prime(p, index);

    const int n_factorial_mod = (n >= static_cast<u64>(p)) ? 0 : factorial_mod_small(n, p);
    return normalize_mod(static_cast<i64>(n_factorial_mod) - static_cast<i64>(A[static_cast<std::size_t>(index)]), p);
}

int crt_combine(const std::vector<int>& residues, const std::vector<int>& moduli) {
    i64 x = 0;
    i64 current_modulus = 1;

    for (std::size_t i = 0; i < residues.size(); ++i) {
        const int m = moduli[i];
        const int a = residues[i];

        const int inv = inverse_mod_prime(static_cast<int>(current_modulus % m), m);
        const int t = static_cast<int>((static_cast<i64>(normalize_mod(static_cast<i64>(a) - x, m)) * inv) % m);
        x += current_modulus * static_cast<i64>(t);
        current_modulus *= static_cast<i64>(m);
    }

    return normalize_mod(x, kModulus);
}

int solve_mod(const u64 n, const bool allow_multithreading, const unsigned requested_threads) {
    const std::vector<int> moduli(std::begin(kPrimeFactors), std::end(kPrimeFactors));
    std::vector<int> residues(moduli.size(), 0);

    const unsigned thread_count = choose_thread_count(allow_multithreading, requested_threads, moduli.size());
    std::atomic<std::size_t> next_idx(0ULL);

    auto worker = [&]() {
        while (true) {
            const std::size_t idx = next_idx.fetch_add(1ULL, std::memory_order_relaxed);
            if (idx >= moduli.size()) {
                break;
            }
            residues[idx] = solve_mod_prime(moduli[idx], n);
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(thread_count);
    for (unsigned t = 0U; t < thread_count; ++t) {
        pool.emplace_back(worker);
    }
    for (std::thread& th : pool) {
        th.join();
    }

    return crt_combine(residues, moduli);
}

struct ExactAB {
    std::vector<i64> A;
    std::vector<i64> B;
    std::vector<i64> fact;
};

ExactAB compute_exact_ab(const int max_n) {
    ExactAB result;
    result.A.assign(static_cast<std::size_t>(max_n) + 1ULL, 0);
    result.B.assign(static_cast<std::size_t>(max_n) + 1ULL, 0);
    result.fact.assign(static_cast<std::size_t>(max_n) + 1ULL, 1);

    std::vector<std::vector<i64>> binom(
        static_cast<std::size_t>(max_n) + 1ULL,
        std::vector<i64>(static_cast<std::size_t>(max_n) + 1ULL, 0));

    for (int n = 0; n <= max_n; ++n) {
        binom[static_cast<std::size_t>(n)][0] = 1;
        binom[static_cast<std::size_t>(n)][static_cast<std::size_t>(n)] = 1;
        for (int k = 1; k < n; ++k) {
            binom[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)] =
                binom[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(k - 1)] +
                binom[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(k)];
        }
    }

    for (int n = 1; n <= max_n; ++n) {
        result.fact[static_cast<std::size_t>(n)] = result.fact[static_cast<std::size_t>(n - 1)] * static_cast<i64>(n);
    }

    for (int n = 0; n <= max_n; ++n) {
        i64 a_n = result.fact[static_cast<std::size_t>(n)];
        for (int k = 0; k < n; ++k) {
            a_n += binom[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)] * result.A[static_cast<std::size_t>(k)];
        }

        i64 s_n = 0;
        for (int k = 0; k <= n; ++k) {
            s_n += result.fact[static_cast<std::size_t>(n)] / result.fact[static_cast<std::size_t>(k)];
        }

        i64 b_n = -s_n;
        for (int k = 0; k < n; ++k) {
            b_n += binom[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)] * result.B[static_cast<std::size_t>(k)];
        }

        result.A[static_cast<std::size_t>(n)] = a_n;
        result.B[static_cast<std::size_t>(n)] = b_n;
    }

    return result;
}

bool run_checkpoints(const Options& options) {
    const ExactAB exact = compute_exact_ab(12);

    if (exact.A[10] != 328161643LL || exact.B[10] != -652694486LL) {
        std::cerr << "Checkpoint failed: a(10) decomposition mismatch. got A(10)="
                  << exact.A[10] << ", B(10)=" << exact.B[10] << '\n';
        return false;
    }

    for (int n = 0; n <= 12; ++n) {
        const i64 expected_sum = exact.A[static_cast<std::size_t>(n)] + exact.B[static_cast<std::size_t>(n)];
        const int expected_mod = normalize_mod(expected_sum, kModulus);
        const int got_mod = solve_mod(static_cast<u64>(n), false, 1U);
        if (got_mod != expected_mod) {
            std::cerr << "Checkpoint failed at n=" << n << ": expected " << expected_mod
                      << ", got " << got_mod << '\n';
            return false;
        }
    }

    for (const int p : kPrimeFactors) {
        const u64 period = static_cast<u64>(p) * static_cast<u64>(p - 1);
        constexpr u64 kSpotChecks = 80ULL;
        const u64 max_index = static_cast<u64>(p) + period + kSpotChecks;
        const std::vector<int> A = compute_A_mod_prime(p, max_index);

        for (u64 t = 0ULL; t < kSpotChecks; ++t) {
            const u64 lhs_idx = static_cast<u64>(p) + t;
            const u64 rhs_idx = lhs_idx + period;
            if (A[static_cast<std::size_t>(lhs_idx)] != A[static_cast<std::size_t>(rhs_idx)]) {
                std::cerr << "Checkpoint failed: period spot-check mismatch for p=" << p
                          << " at offset " << t << '\n';
                return false;
            }
        }
    }

    const unsigned threads = choose_thread_count(options.allow_multithreading,
                                                 options.requested_threads,
                                                 std::size_t{5});
    std::cout << "Checkpoints passed (threads=" << threads << ").\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const int answer = solve_mod(options.n, options.allow_multithreading, options.requested_threads);
    std::cout << "Answer: " << answer << '\n';
    return 0;
}
