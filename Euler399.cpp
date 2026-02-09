#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace {

using u64 = std::uint64_t;
using i64 = long long;
using boost::multiprecision::cpp_dec_float_50;

constexpr u64 kLast16Mod = 10000000000000000ULL;

struct Options {
    i64 target = 100000000;
    i64 k_upper = 220000000;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoll(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_i64_after_prefix(arg, "--target=", options.target) ||
            parse_i64_after_prefix(arg, "--k-upper=", options.k_upper)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 1 && options.k_upper >= 1000;
}

u64 mod_pow_u64(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0) {
        if (e & 1ULL) {
            result = static_cast<u64>((__uint128_t)result * cur % mod);
        }
        cur = static_cast<u64>((__uint128_t)cur * cur % mod);
        e >>= 1ULL;
    }
    return result;
}

std::pair<u64, u64> fib_pair_mod(u64 n, u64 mod) {
    if (n == 0) {
        return {0ULL, 1ULL};
    }
    const auto [a, b] = fib_pair_mod(n >> 1ULL, mod);
    const u64 two_b = static_cast<u64>((2ULL * (__uint128_t)b) % mod);
    const u64 two_b_minus_a = (two_b + mod - a) % mod;
    const u64 c = static_cast<u64>((__uint128_t)a * two_b_minus_a % mod);
    const u64 d = static_cast<u64>(((__uint128_t)a * a + (__uint128_t)b * b) % mod);
    if (n & 1ULL) {
        return {d, (c + d) % mod};
    }
    return {c, d};
}

u64 fibonacci_mod(u64 n, u64 mod) {
    return fib_pair_mod(n, mod).first;
}

int legendre_5_mod_p(const int p) {
    if (p == 5) {
        return 0;
    }
    const u64 ls = mod_pow_u64(5ULL, static_cast<u64>((p - 1) / 2), static_cast<u64>(p));
    if (ls == 1ULL) {
        return 1;
    }
    if (ls == static_cast<u64>(p - 1)) {
        return -1;
    }
    return 0;
}

std::vector<int> build_spf(const int limit, std::vector<int>& primes) {
    std::vector<int> spf(static_cast<std::size_t>(limit + 1), 0);
    primes.clear();
    primes.reserve(limit / 10);
    for (int i = 2; i <= limit; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (const int p : primes) {
            const i64 x = 1LL * i * p;
            if (x > limit) {
                break;
            }
            spf[static_cast<std::size_t>(x)] = p;
            if (p == spf[static_cast<std::size_t>(i)]) {
                break;
            }
        }
    }
    return spf;
}

std::vector<int> unique_prime_factors(int x, const std::vector<int>& spf) {
    std::vector<int> factors;
    int n = x;
    while (n > 1) {
        const int p = spf[static_cast<std::size_t>(n)];
        factors.push_back(p);
        while (n % p == 0) {
            n /= p;
        }
    }
    return factors;
}

int rank_of_apparition_prime(const int p, const std::vector<int>& spf) {
    if (p == 2) {
        return 3;
    }
    if (p == 5) {
        return 5;
    }

    const int leg = legendre_5_mod_p(p);
    int r = p - leg;
    const std::vector<int> factors = unique_prime_factors(r, spf);
    for (const int q : factors) {
        while (r % q == 0) {
            const int candidate = r / q;
            if (fibonacci_mod(static_cast<u64>(candidate), static_cast<u64>(p)) == 0ULL) {
                r = candidate;
            } else {
                break;
            }
        }
    }
    return r;
}

std::vector<int> reduce_periods(std::vector<int> periods) {
    std::sort(periods.begin(), periods.end());
    periods.erase(std::unique(periods.begin(), periods.end()), periods.end());

    std::vector<int> reduced;
    reduced.reserve(periods.size());
    for (const int b : periods) {
        bool redundant = false;
        for (const int kept : reduced) {
            if (kept > b) {
                break;
            }
            if (b % kept == 0) {
                redundant = true;
                break;
            }
        }
        if (!redundant) {
            reduced.push_back(b);
        }
    }
    return reduced;
}

std::vector<int> periods_up_to_k(const i64 k_upper) {
    const int p_limit = static_cast<int>(k_upper / 3);
    std::vector<int> primes;
    const std::vector<int> spf = build_spf(p_limit + 2, primes);

    const unsigned hw = std::thread::hardware_concurrency();
    const int threads = std::max(1U, hw == 0 ? 4U : hw);
    std::vector<std::vector<int>> locals(static_cast<std::size_t>(threads));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));

    auto worker = [&](const int tid) {
        const std::size_t total = primes.size();
        const std::size_t l = total * static_cast<std::size_t>(tid) / static_cast<std::size_t>(threads);
        const std::size_t r = total * static_cast<std::size_t>(tid + 1) / static_cast<std::size_t>(threads);

        std::vector<int>& out = locals[static_cast<std::size_t>(tid)];
        out.reserve((r - l) / 20 + 16);

        for (std::size_t i = l; i < r; ++i) {
            const int p = primes[i];
            const int z = rank_of_apparition_prime(p, spf);
            const i64 period = static_cast<i64>(p) * static_cast<i64>(z);
            if (period <= k_upper) {
                out.push_back(static_cast<int>(period));
            }
        }
    };

    for (int t = 0; t < threads; ++t) {
        workers.emplace_back(worker, t);
    }
    for (std::thread& th : workers) {
        th.join();
    }

    std::vector<int> periods;
    for (const auto& v : locals) {
        periods.insert(periods.end(), v.begin(), v.end());
    }
    return reduce_periods(std::move(periods));
}

i64 nth_squarefree_fib_index(const i64 target, const i64 k_upper, const std::vector<int>& periods) {
    std::vector<std::uint8_t> bad(static_cast<std::size_t>(k_upper + 1), 0U);
    for (const int b : periods) {
        for (i64 x = b; x <= k_upper; x += b) {
            bad[static_cast<std::size_t>(x)] = 1U;
        }
    }

    i64 count = 0;
    for (i64 i = 1; i <= k_upper; ++i) {
        if (bad[static_cast<std::size_t>(i)] == 0U) {
            ++count;
            if (count == target) {
                return i;
            }
        }
    }
    return -1;
}

std::string last16_and_scientific(const i64 index) {
    const u64 last16 = fibonacci_mod(static_cast<u64>(index), kLast16Mod);

    const cpp_dec_float_50 sqrt5 = sqrt(cpp_dec_float_50(5));
    const cpp_dec_float_50 phi = (cpp_dec_float_50(1) + sqrt5) / cpp_dec_float_50(2);
    const cpp_dec_float_50 lg = cpp_dec_float_50(index) * log10(phi) - log10(sqrt5);

    long long exponent = floor(lg).convert_to<long long>();
    cpp_dec_float_50 frac = lg - cpp_dec_float_50(exponent);
    cpp_dec_float_50 mant = pow(cpp_dec_float_50(10), frac);
    cpp_dec_float_50 rounded = floor(mant * 10 + cpp_dec_float_50("0.5")) / 10;
    if (rounded >= 10) {
        rounded /= 10;
        ++exponent;
    }

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(16) << last16 << ','
        << std::fixed << std::setprecision(1) << rounded << 'e' << exponent;
    return oss.str();
}

bool run_checkpoints() {
    const i64 sample_target = 200;
    const i64 sample_upper = 2000;
    const std::vector<int> sample_periods = periods_up_to_k(sample_upper);
    const i64 idx = nth_squarefree_fib_index(sample_target, sample_upper, sample_periods);
    if (idx <= 0) {
        std::cerr << "Checkpoint failed: could not locate 200th index\n";
        return false;
    }

    const std::string sample = last16_and_scientific(idx);
    if (sample != "1608739584170445,9.7e53") {
        std::cerr << "Checkpoint failed for 200th sample, got " << sample << '\n';
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

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::vector<int> periods = periods_up_to_k(options.k_upper);
    i64 index = nth_squarefree_fib_index(options.target, options.k_upper, periods);

    if (index < 0) {
        std::cerr << "Upper bound too small, no index found up to " << options.k_upper << '\n';
        return 3;
    }

    std::cout << last16_and_scientific(index) << '\n';
    return 0;
}
