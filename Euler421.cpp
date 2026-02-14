#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i64 = std::int64_t;

std::string to_string_u128(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v > 0) {
        s.push_back(static_cast<char>('0' + static_cast<int>(v % 10)));
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::vector<int> sieve_primes(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; static_cast<i64>(p) * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) continue;
        for (int q = p * p; q <= n; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    u64 x = a % mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = static_cast<u64>((static_cast<u128>(r) * x) % mod);
        }
        x = static_cast<u64>((static_cast<u128>(x) * x) % mod);
        e >>= 1ULL;
    }
    return r;
}

i64 ext_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

u64 mod_inv(u64 a, u64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = ext_gcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    if (g != 1) return 0;
    i64 r = x % static_cast<i64>(mod);
    if (r < 0) r += static_cast<i64>(mod);
    return static_cast<u64>(r);
}

std::vector<u64> distinct_prime_factors(u64 n, const std::vector<int>& small_primes) {
    std::vector<u64> f;
    u64 x = n;
    for (const int p : small_primes) {
        const u64 pp = static_cast<u64>(p);
        if (pp * pp > x) break;
        if (x % pp != 0) continue;
        f.push_back(pp);
        while (x % pp == 0) x /= pp;
    }
    if (x > 1) f.push_back(x);
    return f;
}

u64 primitive_root_prime(const u64 p, const std::vector<u64>& factors) {
    const u64 phi = p - 1;
    for (u64 g = 2; g < p; ++g) {
        bool ok = true;
        for (const u64 q : factors) {
            if (mod_pow(g, phi / q, p) == 1ULL) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
    return 0;
}

u64 count_roots_leq_rem(const u64 p, const u64 rem, const int d, const std::vector<int>& small_primes) {
    if (d == 1) {
        return (rem == p - 1) ? 1ULL : 0ULL;
    }

    const u64 q = p - 1;
    const u64 mod = q / static_cast<u64>(d);
    const u64 a = 15ULL / static_cast<u64>(d);
    const u64 b = (q / 2ULL) / static_cast<u64>(d);
    const u64 inv_a = mod_inv(a % mod, mod);
    const u64 k0 = static_cast<u64>((static_cast<u128>(b) * inv_a) % mod);

    const auto factors = distinct_prime_factors(q, small_primes);
    const u64 g = primitive_root_prime(p, factors);
    const u64 r0 = mod_pow(g, k0, p);
    const u64 step = mod_pow(g, mod, p);

    u64 count = 0;
    u64 r = r0;
    for (int i = 0; i < d; ++i) {
        if (r <= rem) ++count;
        r = static_cast<u64>((static_cast<u128>(r) * step) % p);
    }
    return count;
}

u128 fast_sum(const u64 n_limit, const int p_limit) {
    const auto primes = sieve_primes(p_limit);
    const auto small_primes = sieve_primes(static_cast<int>(std::sqrt(static_cast<long double>(p_limit)) + 1.0L));

    u128 total = 0;
    for (const int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        if (p == 2ULL) {
            const u64 cnt = n_limit / 2ULL + (n_limit % 2ULL);
            total += static_cast<u128>(2ULL) * cnt;
            continue;
        }

        const u64 qn = n_limit / p;
        const u64 rem = n_limit - qn * p;
        const int d = std::gcd(15, static_cast<int>(p - 1ULL));
        const u64 t = count_roots_leq_rem(p, rem, d, small_primes);
        const u128 cnt = static_cast<u128>(d) * static_cast<u128>(qn) + static_cast<u128>(t);
        total += static_cast<u128>(p) * cnt;
    }
    return total;
}

u64 s_single_bruteforce(const u64 n, const int m) {
    const auto primes = sieve_primes(m);
    u64 sum = 0;
    for (const int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        const u64 r = mod_pow(n % p, 15ULL, p);
        if ((r + 1ULL) % p == 0ULL) sum += p;
    }
    return sum;
}

u128 brute_small(const u64 n_limit, const int m) {
    u128 s = 0;
    for (u64 n = 1; n <= n_limit; ++n) {
        s += s_single_bruteforce(n, m);
    }
    return s;
}

bool run_checkpoints() {
    if (s_single_bruteforce(2ULL, 10) != 3ULL) return false;
    if (s_single_bruteforce(2ULL, 1000) != 345ULL) return false;
    if (s_single_bruteforce(10ULL, 100) != 31ULL) return false;
    if (s_single_bruteforce(10ULL, 1000) != 483ULL) return false;

    if (fast_sum(200ULL, 1000) != brute_small(200ULL, 1000)) return false;
    if (fast_sum(1000ULL, 5000) != brute_small(1000ULL, 5000)) return false;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        std::cerr << "Checkpoint failed\n";
        return 2;
    }

    const u128 ans = fast_sum(100000000000ULL, 100000000);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
