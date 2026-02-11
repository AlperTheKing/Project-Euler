#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;
using u64 = std::uint64_t;

constexpr i64 MOD = 1'000'000'007LL;

i64 mod_pow(i64 base, u64 exp) {
    i64 result = 1 % MOD;
    i64 cur = base % MOD;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = static_cast<i64>((static_cast<i128>(result) * cur) % MOD);
        }
        cur = static_cast<i64>((static_cast<i128>(cur) * cur) % MOD);
        exp >>= 1ULL;
    }
    return result;
}

i64 mod_inv(i64 x) {
    assert(x % MOD != 0);
    return mod_pow((x % MOD + MOD) % MOD, static_cast<u64>(MOD - 2));
}

std::pair<i64, i64> pow_sum(i64 ratio, u64 n) {
    if (n == 0ULL) {
        return {1LL, 0LL};
    }
    if (n == 1ULL) {
        return {ratio % MOD, 1LL};
    }
    if ((n & 1ULL) == 0ULL) {
        auto [p, s] = pow_sum(ratio, n >> 1ULL);
        const i64 p2 = static_cast<i64>((static_cast<i128>(p) * p) % MOD);
        const i64 s2 = static_cast<i64>((static_cast<i128>(s) * ((1LL + p) % MOD)) % MOD);
        return {p2, s2};
    }
    auto [p, s] = pow_sum(ratio, n - 1ULL);
    const i64 p2 = static_cast<i64>((static_cast<i128>(p) * (ratio % MOD)) % MOD);
    const i64 s2 = (s + p) % MOD;
    return {p2, s2};
}

std::pair<i64, u64> quotient_mod_from_exponents(u64 e_prev, u64 e_cur) {
    const u64 a = e_prev / e_cur;
    const u64 b = e_prev % e_cur;
    const i64 two_b = mod_pow(2LL, b);
    const i64 ratio = mod_pow(2LL, e_cur);
    const i64 geom = pow_sum(ratio, a).second;
    const i64 q_mod = static_cast<i64>((static_cast<i128>(two_b) * geom) % MOD);
    return {q_mod, b};
}

i64 inverse_mersenne_mod_mod(u64 u, u64 w) {
    u64 e_prev = u;
    u64 e_cur = w;
    i64 t_prev = 0;
    i64 t_cur = 1;
    int steps = 0;

    while (e_cur != 0ULL) {
        auto [q_mod, b] = quotient_mod_from_exponents(e_prev, e_cur);
        i64 t_next = static_cast<i64>((t_prev - static_cast<i128>(q_mod) * t_cur) % MOD);
        if (t_next < 0) {
            t_next += MOD;
        }
        e_prev = e_cur;
        e_cur = b;
        t_prev = t_cur;
        t_cur = t_next;
        ++steps;
    }

    const i64 a_mod = (mod_pow(2LL, u) - 1 + MOD) % MOD;
    if (steps % 2 == 1) {
        return t_prev;
    }
    return (a_mod + t_prev) % MOD;
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

u64 inverse_u64_mod(u64 a, u64 m) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = ext_gcd(static_cast<i64>(a), static_cast<i64>(m), x, y);
    assert(g == 1);
    i64 r = x % static_cast<i64>(m);
    if (r < 0) {
        r += static_cast<i64>(m);
    }
    return static_cast<u64>(r);
}

u64 pow_u64(u64 base, int exp) {
    u64 result = 1ULL;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

u64 p_general(u64 a, u64 b) {
    assert(a <= b);
    assert(std::gcd(a, b) == 1ULL);

    const u64 inv = inverse_u64_mod(a % b, b);
    const u64 y = static_cast<u64>((static_cast<i128>(a) * inv - 1) / static_cast<i128>(b));
    const u64 s1 = inv + y - 1;
    const u64 s2 = (b - inv) + (a - y) - 1;
    return 2ULL * std::min(s1, s2);
}

i64 p_mersenne_from_exponents(u64 u, u64 v) {
    assert(u < v);
    assert(std::gcd(u, v) == 1ULL);

    const i64 a_mod = (mod_pow(2LL, u) - 1 + MOD) % MOD;
    const i64 b_mod = (mod_pow(2LL, v) - 1 + MOD) % MOD;
    const i64 inv_a_mod = mod_inv(a_mod);

    const u64 w = v % u;
    const i64 r_mod = inverse_mersenne_mod_mod(u, w);
    const u64 t = inverse_u64_mod(w, u);

    i64 d_mod = 0;
    i64 x_mod = 0;
    if (2ULL * t <= u) {
        d_mod = r_mod;
        i64 rhs = (static_cast<i64>((static_cast<i128>(b_mod) * d_mod) % MOD) - 1 + MOD) % MOD;
        x_mod = static_cast<i64>((static_cast<i128>(rhs) * inv_a_mod) % MOD);
    } else {
        d_mod = (a_mod - r_mod + MOD) % MOD;
        i64 rhs = (static_cast<i64>((static_cast<i128>(b_mod) * d_mod) % MOD) + 1) % MOD;
        x_mod = static_cast<i64>((static_cast<i128>(rhs) * inv_a_mod) % MOD);
    }

    const i64 m_mod = (x_mod + d_mod - 1 + MOD) % MOD;
    return static_cast<i64>((2LL * m_mod) % MOD);
}

std::vector<int> primes_below(int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n), true);
    is_prime[0] = false;
    is_prime[1] = false;
    for (int p = 2; static_cast<i64>(p) * p < n; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        for (int x = p * p; x < n; x += p) {
            is_prime[x] = false;
        }
    }
    std::vector<int> primes;
    for (int x = 2; x < n; ++x) {
        if (is_prime[x]) {
            primes.push_back(x);
        }
    }
    return primes;
}

}  // namespace

int main() {
    assert(p_general(3, 5) == 4ULL);
    assert(p_general(7, 31) == 20ULL);
    assert(p_general(1234, 4321) == 2780ULL);

    assert(p_mersenne_from_exponents(2, 5) == 20);
    assert(p_mersenne_from_exponents(3, 5) == 20);
    assert(p_mersenne_from_exponents(5, 8) == 164);
    for (u64 u = 2; u <= 9; ++u) {
        for (u64 v = u + 1; v <= 10; ++v) {
            if (std::gcd(u, v) != 1ULL) {
                continue;
            }
            const u64 a = (1ULL << u) - 1ULL;
            const u64 b = (1ULL << v) - 1ULL;
            assert(p_mersenne_from_exponents(u, v) == static_cast<i64>(p_general(a, b)));
        }
    }

    const std::vector<int> primes = primes_below(1000);
    std::vector<u64> p5;
    p5.reserve(primes.size());
    for (int p : primes) {
        p5.push_back(pow_u64(static_cast<u64>(p), 5));
    }

    i64 answer = 0;
    for (std::size_t i = 0; i < p5.size(); ++i) {
        for (std::size_t j = i + 1; j < p5.size(); ++j) {
            answer += p_mersenne_from_exponents(p5[i], p5[j]);
            if (answer >= MOD) {
                answer -= MOD;
            }
        }
    }

    std::cout << answer << '\n';
    return 0;
}
