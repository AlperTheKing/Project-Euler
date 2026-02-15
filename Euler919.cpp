#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

struct Pair {
    i64 k;
    i64 l;
};

static inline i64 isqrt(u64 x) {
    const i64 approx = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    i64 r = approx;
    while ((r + 1) * (r + 1) <= static_cast<u64>(x)) {
        ++r;
    }
    while (r * r > static_cast<u64>(x)) {
        --r;
    }
    return r;
}

i64 fortunate_sum_fast(i64 perimeter_limit) {
    constexpr std::array<Pair, 8> coeffs{{
        {1, 15},
        {3, 5},
        {5, 3},
        {15, 1},
        {2, 30},
        {6, 10},
        {10, 6},
        {30, 2},
    }};

    i64 ans = 0;

    for (const auto [k, l] : coeffs) {
        for (i64 r = 1;; ++r) {
            const i64 sm = k * r * r;
            if (sm > 4LL * perimeter_limit) {
                break;
            }

            for (i64 s = 1;; ++s) {
                const i64 df = l * s * s;
                if (df >= sm) {
                    break;
                }
                if (std::gcd(r, s) != 1) {
                    continue;
                }
                if ((sm + df) & 1LL) {
                    continue;
                }

                const i64 x = (sm + df) / 2;
                const i64 y = (sm - df) / 2;
                const i64 z2 = (x * x - y * y) / 15;

                const i64 z = isqrt(static_cast<u64>(z2));
                if (z * z != z2) {
                    continue;
                }
                if (std::gcd(x, std::gcd(y, z)) != 1) {
                    continue;
                }

                i64 a = x;
                i64 b = y + z;
                i64 c = z;
                while ((a & 3LL) || (b & 3LL)) {
                    a <<= 1;
                    b <<= 1;
                    c <<= 1;
                }
                a >>= 2;
                b >>= 2;
                if (std::gcd(a, std::gcd(b, c)) <= 1) {
                    if (a + b + c > 2LL * std::max(a, std::max(b, c)) && b >= c) {
                        const i64 p = a + b + c;
                        const i64 ct = perimeter_limit / p;
                        ans += p * ct * (ct + 1) / 2;
                    }
                }

                if (z >= y) {
                    continue;
                }

                a = x;
                b = y - z;
                c = z;
                while ((a & 3LL) || (b & 3LL)) {
                    a <<= 1;
                    b <<= 1;
                    c <<= 1;
                }
                a >>= 2;
                b >>= 2;

                if (std::gcd(a, std::gcd(b, c)) <= 1) {
                    if (a + b + c > 2LL * std::max(a, std::max(b, c)) && b >= c) {
                        const i64 p = a + b + c;
                        const i64 ct = perimeter_limit / p;
                        ans += p * ct * (ct + 1) / 2;
                    }
                }
            }
        }
    }

    return ans;
}

}  // namespace

int main() {
    std::cout << fortunate_sum_fast(10'000'000) << '\n';
    return 0;
}
