#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using boost::multiprecision::cpp_int;

constexpr u64 kMod = 999'999'001ULL;

u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mod_mul(result, cur, mod);
        }
        cur = mod_mul(cur, cur, mod);
        exp >>= 1ULL;
    }
    return result;
}

u64 mod_inv(u64 x, u64 mod) {
    return mod_pow(x, mod - 2ULL, mod);
}

std::vector<int> build_spf(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 5);

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > static_cast<u64>(n) || p > spf[i]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }
    return spf;
}

std::vector<std::pair<int, u64>> prime_exponents_superduper(int n) {
    const std::vector<int> spf = build_spf(n);
    std::vector<u64> exp(static_cast<std::size_t>(n + 1), 0ULL);

    for (int t = 1; t <= n; ++t) {
        const u64 coef = static_cast<u64>(n - t + 1) * static_cast<u64>(n - t + 2) / 2ULL;
        int x = t;
        while (x > 1) {
            const int p = spf[static_cast<std::size_t>(x)];
            int cnt = 0;
            while (x % p == 0) {
                x /= p;
                ++cnt;
            }
            exp[static_cast<std::size_t>(p)] += coef * static_cast<u64>(cnt);
        }
    }

    std::vector<std::pair<int, u64>> factors;
    for (int p = 2; p <= n; ++p) {
        if (exp[static_cast<std::size_t>(p)] != 0ULL) {
            factors.push_back({p, exp[static_cast<std::size_t>(p)]});
        }
    }
    return factors;
}

u64 D_mod(int n, int m, u64 mod) {
    const auto factors = prime_exponents_superduper(n);
    std::vector<u64> dp(static_cast<std::size_t>(m), 0ULL);
    std::vector<u64> next(static_cast<std::size_t>(m), 0ULL);
    std::vector<u64> coeff(static_cast<std::size_t>(m), 0ULL);
    dp[0] = 1ULL;

    for (const auto& [p_int, e] : factors) {
        const u64 p = static_cast<u64>(p_int);
        const u64 q = mod_pow(p, static_cast<u64>(m), mod);

        u64 p_pow_t = 1ULL;
        for (int t = 0; t < m; ++t) {
            if (e >= static_cast<u64>(t)) {
                const u64 cnt = (e - static_cast<u64>(t)) / static_cast<u64>(m) + 1ULL;
                u64 geom = 0ULL;
                if (q == 1ULL) {
                    geom = cnt % mod;
                } else {
                    const u64 num = (mod_pow(q, cnt, mod) + mod - 1ULL) % mod;
                    const u64 den_inv = mod_inv((q + mod - 1ULL) % mod, mod);
                    geom = mod_mul(num, den_inv, mod);
                }
                coeff[static_cast<std::size_t>(t)] = mod_mul(p_pow_t, geom, mod);
            } else {
                coeff[static_cast<std::size_t>(t)] = 0ULL;
            }
            p_pow_t = mod_mul(p_pow_t, p, mod);
        }

        std::fill(next.begin(), next.end(), 0ULL);
        for (int r = 0; r < m; ++r) {
            if (dp[static_cast<std::size_t>(r)] == 0ULL) {
                continue;
            }
            const u64 base = dp[static_cast<std::size_t>(r)];
            for (int t = 0; t < m; ++t) {
                const int idx = (r + t >= m) ? (r + t - m) : (r + t);
                next[static_cast<std::size_t>(idx)] =
                    (next[static_cast<std::size_t>(idx)] + mod_mul(base, coeff[static_cast<std::size_t>(t)], mod)) % mod;
            }
        }
        dp.swap(next);
    }

    return dp[0];
}

cpp_int D_exact_small(int n, int m) {
    const auto factors = prime_exponents_superduper(n);
    std::vector<cpp_int> dp(static_cast<std::size_t>(m), cpp_int(0));
    std::vector<cpp_int> next(static_cast<std::size_t>(m), cpp_int(0));
    std::vector<cpp_int> coeff(static_cast<std::size_t>(m), cpp_int(0));
    dp[0] = 1;

    for (const auto& [p_int, e] : factors) {
        const cpp_int p = p_int;

        for (int t = 0; t < m; ++t) {
            if (e < static_cast<u64>(t)) {
                coeff[static_cast<std::size_t>(t)] = 0;
                continue;
            }
            cpp_int term = 1;
            for (int i = 0; i < t; ++i) {
                term *= p;
            }
            cpp_int sum = 0;
            for (u64 a = static_cast<u64>(t); a <= e; a += static_cast<u64>(m)) {
                sum += term;
                for (int i = 0; i < m; ++i) {
                    term *= p;
                }
                if (e - a < static_cast<u64>(m)) {
                    break;
                }
            }
            coeff[static_cast<std::size_t>(t)] = sum;
        }

        std::fill(next.begin(), next.end(), cpp_int(0));
        for (int r = 0; r < m; ++r) {
            if (dp[static_cast<std::size_t>(r)] == 0) {
                continue;
            }
            for (int t = 0; t < m; ++t) {
                const int idx = (r + t >= m) ? (r + t - m) : (r + t);
                next[static_cast<std::size_t>(idx)] += dp[static_cast<std::size_t>(r)] * coeff[static_cast<std::size_t>(t)];
            }
        }
        dp.swap(next);
    }

    return dp[0];
}

void run_validations() {
    const cpp_int sample_exact = D_exact_small(6, 6);
    assert(sample_exact == cpp_int("6368195719791280"));
    const u64 sample_mod = D_mod(6, 6, kMod);
    assert(sample_mod == static_cast<u64>(sample_exact % kMod));
}

}  // namespace

int main() {
    run_validations();
    std::cout << D_mod(1'000, 1'000, kMod) << '\n';
    return 0;
}
