#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'000ULL;

u64 isqrt_u64(const u64 n) {
    const long double approx = static_cast<long double>(n);
    u64 x = static_cast<u64>(std::sqrt(approx));
    while ((x + 1ULL) * (x + 1ULL) <= n) {
        ++x;
    }
    while (x * x > n) {
        --x;
    }
    return x;
}

u64 max_r(const u64 n_limit) {
    const u128 rhs = 2ULL * static_cast<u128>(n_limit) * static_cast<u128>(n_limit);
    u64 lo = 0ULL;
    u64 hi = static_cast<u64>(std::sqrt(static_cast<long double>(n_limit))) + 10ULL;

    while (13ULL * static_cast<u128>(hi) * hi * hi * hi <= rhs) {
        hi *= 2ULL;
    }

    while (lo + 1ULL < hi) {
        const u64 mid = lo + (hi - lo) / 2ULL;
        const u128 lhs = 13ULL * static_cast<u128>(mid) * mid * mid * mid;
        if (lhs <= rhs) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return lo;
}

u128 S(const u64 n_limit) {
    const u64 r_limit = max_r(n_limit);
    const u64 m_limit = isqrt_u64(r_limit);

    u128 total = 0;

    for (u64 m = 1ULL; m <= m_limit; ++m) {
        for (u64 n = 0ULL; n < m; ++n) {
            if (((m - n) & 1ULL) == 0ULL) {
                continue;
            }
            if (std::gcd(m, n) != 1ULL) {
                continue;
            }

            const u64 r = m * m + n * n;
            if (r > r_limit) {
                continue;
            }

            const std::int64_t u = static_cast<std::int64_t>(m * m - n * n);
            const std::int64_t v = static_cast<std::int64_t>(2ULL * m * n);

            std::pair<u64, u64> cand[2] = {
                {static_cast<u64>(std::llabs(3LL * u - 2LL * v)), static_cast<u64>(std::llabs(2LL * u + 3LL * v))},
                {static_cast<u64>(std::llabs(3LL * u + 2LL * v)), static_cast<u64>(std::llabs(2LL * u - 3LL * v))},
            };

            for (int i = 0; i < 2; ++i) {
                if (i == 1 && cand[1] == cand[0]) {
                    continue;
                }

                u64 p = cand[i].first;
                u64 q = cand[i].second;
                if (p == 0ULL || q == 0ULL) {
                    continue;
                }
                if (p < q) {
                    std::swap(p, q);
                }
                if (std::gcd(p, q) != 1ULL) {
                    continue;
                }

                const u128 x = static_cast<u128>(q) * static_cast<u128>(r);
                const u128 y = static_cast<u128>(p) * static_cast<u128>(r);
                const u128 z = static_cast<u128>(p) * static_cast<u128>(q);
                if (x > n_limit || y > n_limit || z > n_limit) {
                    continue;
                }

                total += x + y + z;
            }
        }
    }

    return total;
}

u64 to_u64(const u128 value) {
    return static_cast<u64>(value);
}

}  // namespace

int main() {
    assert(to_u64(S(100ULL)) == 124ULL);
    assert(to_u64(S(1'000ULL)) == 1'470ULL);
    assert(to_u64(S(100'000ULL)) == 2'340'084ULL);

    const u64 answer = static_cast<u64>(S(10'000'000'000'000'000ULL) % kMod);
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}
