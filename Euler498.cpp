#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kPrime = 999'999'937ULL;

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        e >>= 1ULL;
    }
    return result;
}

u64 binom_mod_small_prime_digit(u64 n, u64 k, const u64 prime) {
    if (k > n) {
        return 0ULL;
    }
    k = std::min(k, n - k);
    if (k == 0ULL) {
        return 1ULL;
    }

    u64 numerator = 1ULL;
    u64 denominator = 1ULL;
    for (u64 i = 1ULL; i <= k; ++i) {
        numerator = static_cast<u64>((static_cast<u128>(numerator) * (n - k + i)) % prime);
        denominator = static_cast<u64>((static_cast<u128>(denominator) * i) % prime);
    }
    const u64 inv_denominator = mod_pow(denominator, prime - 2ULL, prime);
    return static_cast<u64>((static_cast<u128>(numerator) * inv_denominator) % prime);
}

u64 binom_mod_prime_lucas(u64 n, u64 k, const u64 prime) {
    if (k > n) {
        return 0ULL;
    }
    u64 result = 1ULL;
    u64 nn = n;
    u64 kk = k;
    while (nn > 0ULL || kk > 0ULL) {
        const u64 ni = nn % prime;
        const u64 ki = kk % prime;
        if (ki > ni) {
            return 0ULL;
        }
        const u64 digit = binom_mod_small_prime_digit(ni, ki, prime);
        result = static_cast<u64>((static_cast<u128>(result) * digit) % prime);
        nn /= prime;
        kk /= prime;
    }
    return result;
}

u64 coefficient_mod(const u64 n, const u64 m, const u64 d, const u64 prime) {
    if (d >= m || m > n) {
        return 0ULL;
    }
    const u64 first = binom_mod_prime_lucas(n, d, prime);
    const u64 second = binom_mod_prime_lucas(n - d - 1ULL, m - d - 1ULL, prime);
    return static_cast<u64>((static_cast<u128>(first) * second) % prime);
}

u64 binom_exact_small(const u64 n, const u64 k) {
    if (k > n) {
        return 0ULL;
    }
    const u64 kk = std::min(k, n - k);
    u128 result = 1;
    for (u64 i = 1ULL; i <= kk; ++i) {
        result = (result * static_cast<u128>(n - kk + i)) / static_cast<u128>(i);
    }
    return static_cast<u64>(result);
}

u64 coefficient_abs_exact_small(const u64 n, const u64 m, const u64 d) {
    if (d >= m || m > n) {
        return 0ULL;
    }
    const u64 first = binom_exact_small(n, d);
    const u64 second = binom_exact_small(n - d - 1ULL, m - d - 1ULL);
    return first * second;
}

u64 coefficient_abs_bruteforce_small(const u64 n, const u64 m, const u64 d) {
    if (d >= m || m > n) {
        return 0ULL;
    }

    std::vector<long long> coeffs(static_cast<std::size_t>(m), 0LL);
    for (u64 k = 0ULL; k < m; ++k) {
        const u64 choose_n_k = binom_exact_small(n, k);
        for (u64 j = 0ULL; j <= k; ++j) {
            const long long choose_k_j = static_cast<long long>(binom_exact_small(k, j));
            const long long sign = ((k - j) % 2ULL == 0ULL) ? 1LL : -1LL;
            coeffs[static_cast<std::size_t>(j)] +=
                static_cast<long long>(choose_n_k) * choose_k_j * sign;
        }
    }
    const long long value = coeffs[static_cast<std::size_t>(d)];
    return static_cast<u64>(value >= 0LL ? value : -value);
}

bool run_checkpoints() {
    if (coefficient_abs_exact_small(6ULL, 3ULL, 1ULL) != 24ULL) {
        std::cerr << "Checkpoint failed: C(6,3,1)\n";
        return false;
    }
    if (coefficient_abs_exact_small(100ULL, 10ULL, 4ULL) != 227'197'811'615'775ULL) {
        std::cerr << "Checkpoint failed: C(100,10,4)\n";
        return false;
    }

    for (u64 n = 1ULL; n <= 12ULL; ++n) {
        for (u64 m = 1ULL; m <= n; ++m) {
            for (u64 d = 0ULL; d < m; ++d) {
                const u64 closed_form = coefficient_abs_exact_small(n, m, d);
                const u64 brute = coefficient_abs_bruteforce_small(n, m, d);
                if (closed_form != brute) {
                    std::cerr << "Bruteforce mismatch at n=" << n << ", m=" << m << ", d=" << d
                              << '\n';
                    return false;
                }
                const u64 via_mod = coefficient_mod(n, m, d, kPrime);
                if (via_mod != (closed_form % kPrime)) {
                    std::cerr << "Modulo mismatch at n=" << n << ", m=" << m << ", d=" << d
                              << '\n';
                    return false;
                }
            }
        }
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 10'000'000'000'000ULL;
    constexpr u64 m = 1'000'000'000'000ULL;
    constexpr u64 d = 10'000ULL;

    const u64 answer = coefficient_mod(n, m, d, kPrime);
    std::cout << answer << '\n';
    return 0;
}
