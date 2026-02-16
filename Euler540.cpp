#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static inline u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((u128)(r + 1) * (r + 1) <= (u128)x) ++r;
    while ((u128)r * r > (u128)x) --r;
    return r;
}

static inline u64 icbrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while ((u128)(r + 1) * (r + 1) * (r + 1) <= (u128)x) ++r;
    while ((u128)r * r * r > (u128)x) --r;
    return r;
}

static i64 helper(const u64 L) {
    i64 res = 0;
    u64 m0 = static_cast<u64>(0.5L + std::sqrt(0.25L + static_cast<long double>(L - 1) / 2.0L));
    u64 K = m0 / 2;
    res += static_cast<i64>(K * K);

    u64 m = 2 * K + 1;
    u64 n = isqrt_u64((L > m * m) ? (L - m * m) : 0ULL);
    while (true) {
        u64 z = (L > m * m) ? (L - m * m) : 0ULL;
        if (z == 0ULL) {
            break;
        }
        while ((u128)n * n > (u128)z) {
            --n;
        }
        res += static_cast<i64>(n / 2);
        ++m;

        z = (L > m * m) ? (L - m * m) : 0ULL;
        if (z == 0ULL) {
            break;
        }
        while ((u128)n * n > (u128)z) {
            --n;
        }
        res += static_cast<i64>((n + 1) / 2);
        ++m;
    }

    return res;
}

static u64 solve(const u64 N) {
    const u64 L = icbrt_u64(N);
    std::vector<i64> v(static_cast<std::size_t>(L + 1), 0);
    std::vector<i64> bigV(static_cast<std::size_t>(L / 2 + 1), 0);

    for (u64 x = 1; x <= L; ++x) {
        i64 res = helper(x);
        const u64 c = icbrt_u64(x);

        for (u64 g = 3; g <= c; g += 2) {
            res -= v[static_cast<std::size_t>(x / (g * g))];
        }

        i64 prev_y = static_cast<i64>((isqrt_u64(x) - 1) / 2);
        const u64 z_end = c + ((x / (c * c) != c) ? 1ULL : 0ULL);
        for (u64 z = 1; z < z_end; ++z) {
            i64 y = static_cast<i64>((isqrt_u64(x / (z + 1)) - 1) / 2);
            res -= (prev_y - y) * v[static_cast<std::size_t>(z)];
            prev_y = y;
        }
        v[static_cast<std::size_t>(x)] = res;
    }

    for (i64 a = static_cast<i64>((L - 1) / 2); a >= 0; --a) {
        const u64 x = static_cast<u64>(2 * a + 1);
        const u64 k = N / (x * x);
        i64 res = helper(k);
        const u64 c = icbrt_u64(k);

        for (u64 b = 1; b <= (c - 1) / 2; ++b) {
            const u64 g = 2 * b + 1;
            const u64 k_gg = k / (g * g);
            if (k_gg <= L) {
                res -= v[static_cast<std::size_t>(k_gg)];
            } else {
                const u64 idx = 2ULL * static_cast<u64>(a) * b + static_cast<u64>(a) + b;
                res -= bigV[static_cast<std::size_t>(idx)];
            }
        }

        i64 prev_y = static_cast<i64>((isqrt_u64(k) - 1) / 2);
        const u64 z_end = c + ((k / (c * c) != c) ? 1ULL : 0ULL);
        for (u64 z = 1; z < z_end; ++z) {
            i64 y = static_cast<i64>((isqrt_u64(k / (z + 1)) - 1) / 2);
            res -= (prev_y - y) * v[static_cast<std::size_t>(z)];
            prev_y = y;
        }

        bigV[static_cast<std::size_t>(a)] = res;
    }

    return static_cast<u64>(bigV[0]);
}

int main() {
    assert(solve(20ULL) == 3ULL);
    assert(solve(1'000'000ULL) == 159139ULL);
    std::cout << solve(3'141'592'653'589'793ULL) << '\n';
    return 0;
}
