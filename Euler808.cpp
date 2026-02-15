#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static inline u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod(r, a, mod);
        }
        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    static constexpr u64 WITNESSES[] = {2ULL, 325ULL, 9'375ULL, 28'178ULL, 450'775ULL, 9'780'504ULL, 1'795'265'022ULL};
    for (u64 a : WITNESSES) {
        if (a % n == 0) {
            continue;
        }
        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) {
            continue;
        }
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                comp = false;
                break;
            }
        }
        if (comp) {
            return false;
        }
    }
    return true;
}

static std::vector<u32> sieve_primes(const u32 limit) {
    std::vector<bool> is_prime_vec(static_cast<std::size_t>(limit + 1), true);
    if (limit >= 0) {
        is_prime_vec[0] = false;
    }
    if (limit >= 1) {
        is_prime_vec[1] = false;
    }
    for (u32 p = 2; static_cast<u64>(p) * p <= limit; ++p) {
        if (!is_prime_vec[p]) {
            continue;
        }
        for (u32 q = p * p; q <= limit; q += p) {
            is_prime_vec[q] = false;
        }
    }
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(limit / 10));
    for (u32 i = 2; i <= limit; ++i) {
        if (is_prime_vec[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

static u64 reverse_digits(u64 x) {
    u64 r = 0;
    while (x > 0) {
        r = r * 10 + (x % 10);
        x /= 10;
    }
    return r;
}

static bool is_palindrome(u64 x) {
    return x == reverse_digits(x);
}

static u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

static std::vector<u64> first_reversible_prime_squares(int need) {
    for (u32 limit = 1U << 16U;; limit <<= 1U) {
        const auto primes = sieve_primes(limit);
        std::vector<u64> vals;
        vals.reserve(static_cast<std::size_t>(need + 16));

        for (u32 p : primes) {
            const u64 sq = static_cast<u64>(p) * static_cast<u64>(p);
            if (is_palindrome(sq)) {
                continue;
            }
            const u64 rev = reverse_digits(sq);
            const u64 r = isqrt_u64(rev);
            if (r * r == rev && is_prime(r)) {
                vals.push_back(sq);
                if (static_cast<int>(vals.size()) >= need) {
                    vals.resize(static_cast<std::size_t>(need));
                    return vals;
                }
            }
        }
        if (limit >= (1U << 30U)) {
            break;
        }
    }

    return {};
}

int main() {
    const auto vals = first_reversible_prime_squares(50);

    assert(vals[0] == 169ULL);
    assert(vals[1] == 961ULL);

    u64 sum = 0;
    for (u64 v : vals) {
        sum += v;
    }

    std::cout << sum << '\n';
    return 0;
}
