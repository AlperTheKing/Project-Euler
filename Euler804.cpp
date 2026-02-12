#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

using i64 = std::int64_t;
using u64 = std::uint64_t;

static u64 isqrt_u64(u64 n) {
    long double x = std::sqrt(static_cast<long double>(n));
    u64 r = static_cast<u64>(x);
    while ((r + 1) > r && (r + 1) * (r + 1) <= n) {
        ++r;
    }
    while (r * r > n) {
        --r;
    }
    return r;
}

static i64 solve(u64 N) {
    const u64 y_max = isqrt_u64((4ULL * N) / 163ULL);
    i64 total = 0;

    for (i64 y = -static_cast<i64>(y_max); y <= static_cast<i64>(y_max); ++y) {
        const u64 yy = static_cast<u64>(y >= 0 ? y : -y);
        const u64 rem = 4ULL * N - 163ULL * yy * yy;
        const u64 m = isqrt_u64(rem);

        i64 cnt;
        if ((y & 1LL) == 0) {
            cnt = (m & 1ULL) == 0ULL ? static_cast<i64>(m + 1) : static_cast<i64>(m);
        } else {
            cnt = (m & 1ULL) == 0ULL ? static_cast<i64>(m) : static_cast<i64>(m + 1);
        }

        total += cnt;
    }

    return total - 1;
}

int main() {
    assert(solve(1'000ULL) == 474);
    assert(solve(1'000'000ULL) == 492'128);

    std::cout << solve(10'000'000'000'000'000ULL) << '\n';
    return 0;
}
