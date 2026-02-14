#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr u64 kMod = 1'000'000'007ULL;

static inline i64 isqrt_floor(const i64 n) {
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

static inline i64 isqrt_ceil(const i64 n) {
    const i64 r = isqrt_floor(n);
    return (r * r == n) ? r : (r + 1);
}

static i128 psi_prefix_exact(const i64 n) {
    i128 ans = static_cast<i128>(n + 18) * static_cast<i128>(n - 1) * static_cast<i128>(n - 2) / 6;

    for (i64 x = 1; 4 * x * x <= n; ++x) {
        for (i64 y = x; 4 * x * y <= n; ++y) {
            const i64 disc = 4 * x * y * (x + 1) * (y + 1);
            const i64 z_min = isqrt_ceil(disc) + 2 * x * y;
            const i64 z_max = n - 1 - x - y;
            if (z_min > z_max) break;

            i64 cnt = 2 * (z_max - z_min + 1);
            if (z_min * z_min == 4 * (x + y + z_min + 1) * x * y) --cnt;
            ans += static_cast<i128>(3) * static_cast<i128>(cnt) * static_cast<i128>((x == y) ? 1 : 2);
        }
    }
    return ans;
}

static u64 psi_prefix_mod(const i64 n) {
    return static_cast<u64>(psi_prefix_exact(n) % kMod);
}

int main() {
    assert(static_cast<u64>(psi_prefix_exact(10)) == 345ULL);
    assert(static_cast<u64>(psi_prefix_exact(1000)) == 172'166'601ULL);
    std::cout << psi_prefix_mod(100'000'000) << '\n';
    return 0;
}
