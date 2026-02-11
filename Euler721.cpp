#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

constexpr u64 kMod = 999'999'937ULL;

struct Pair {
    u64 x;
    u64 y;
};

u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(__builtin_sqrtl(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

u64 ceil_sqrt_u64(const u64 n) {
    const u64 r = isqrt_u64(n);
    return (r * r == n) ? r : (r + 1ULL);
}

u64 pow_u64_mod(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL % mod;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = static_cast<u64>(static_cast<i128>(result) * base % mod);
        }
        base = static_cast<u64>(static_cast<i128>(base) * base % mod);
        exp >>= 1ULL;
    }
    return result;
}

Pair mul_pair(const Pair a, const Pair b, const u64 a_mod) {
    const u64 real = static_cast<u64>((static_cast<i128>(a.x) * b.x +
                                       static_cast<i128>(a.y) * b.y % kMod * a_mod) %
                                      kMod);
    const u64 imag = static_cast<u64>((static_cast<i128>(a.x) * b.y + static_cast<i128>(a.y) * b.x) % kMod);
    return {real, imag};
}

u64 lucas_sum_mod(const u64 a, const u64 m, u64 n) {
    const u64 a_mod = a % kMod;
    Pair base{m % kMod, 1ULL};
    Pair result{1ULL, 0ULL};

    while (n > 0ULL) {
        if (n & 1ULL) {
            result = mul_pair(result, base, a_mod);
        }
        n >>= 1ULL;
        if (n > 0ULL) {
            base = mul_pair(base, base, a_mod);
        }
    }

    return (2ULL * result.x) % kMod;
}

u64 f_mod(const u64 a, const u64 n) {
    const u64 m = ceil_sqrt_u64(a);
    const bool is_square = (m * m == a);
    u64 value = lucas_sum_mod(a, m, n);
    if (!is_square) {
        value = (value + kMod - 1ULL) % kMod;
    }
    return value;
}

u128 f_exact_small(const u64 a, const u64 n) {
    const u64 m = ceil_sqrt_u64(a);
    const u64 d = m * m - a;

    if (d == 0ULL) {
        u128 p = 1;
        const u128 base = static_cast<u128>(2ULL * m);
        for (u64 i = 0ULL; i < n; ++i) {
            p *= base;
        }
        return p;
    }

    u128 s0 = 2;
    u128 s1 = static_cast<u128>(2ULL * m);
    if (n == 0ULL) {
        return 1;
    }
    if (n == 1ULL) {
        return s1 - 1;
    }

    for (u64 i = 2ULL; i <= n; ++i) {
        const u128 s = static_cast<u128>(2ULL * m) * s1 - static_cast<u128>(d) * s0;
        s0 = s1;
        s1 = s;
    }
    return s1 - 1;
}

u64 G_mod(const int limit) {
    u64 sum = 0ULL;
    u64 m = 1ULL;
    u64 sq = 1ULL;

    for (u64 a = 1ULL; a <= static_cast<u64>(limit); ++a) {
        while (sq < a) {
            ++m;
            sq = m * m;
        }
        const bool is_square = (sq == a);
        u64 value = lucas_sum_mod(a, m, a * a);
        if (!is_square) {
            value = (value + kMod - 1ULL) % kMod;
        }

        sum += value;
        if (sum >= kMod) {
            sum -= kMod;
        }
    }

    return sum;
}

}  // namespace

int main() {
    assert(f_exact_small(5, 2) == 27);
    assert(f_exact_small(5, 5) == 3935);
    assert(f_mod(5, 2) == 27 % kMod);
    assert(f_mod(5, 5) == 3935 % kMod);
    assert(G_mod(1000) == 163'861'845ULL);

    std::cout << G_mod(5'000'000) << '\n';
    return 0;
}
