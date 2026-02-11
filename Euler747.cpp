#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

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

u128 easy_sum(const u64 m) {
    const u128 mm = static_cast<u128>(m);
    return (mm * mm * mm + 15ULL * mm * mm - 52ULL * mm + 36ULL) / 6ULL;
}

u128 hard_one_corner(const u64 m) {
    u128 total = 0;

    for (u64 x = 1ULL;; ++x) {
        const u64 a = x * (x + 1ULL);
        const u64 nmin_diag = 4ULL * x * x + 4ULL * x + 1ULL;
        if (nmin_diag > m) {
            break;
        }

        u64 y = x;
        u64 d = a * y * (y + 1ULL);
        u64 d4 = 4ULL * d;
        u64 s = a;
        u64 s4 = 2ULL * a;

        for (;; ++y) {
            while ((s4 + 1ULL) * (s4 + 1ULL) <= d4) {
                ++s4;
            }
            while (s4 * s4 > d4) {
                --s4;
            }

            u64 nmin = 2ULL * x * y + x + y + 1ULL + s4;
            if (s4 * s4 < d4) {
                ++nmin;
            }
            if (nmin > m) {
                break;
            }

            while ((s + 1ULL) * (s + 1ULL) <= d) {
                ++s;
            }
            while (s * s > d) {
                --s;
            }

            const bool is_square = (s * s == d);
            const u64 ways = 2ULL * (m - nmin + 1ULL) - static_cast<u64>(is_square);
            total += (x == y) ? static_cast<u128>(ways) : static_cast<u128>(2ULL * ways);

            const u64 y_next = y + 1ULL;
            const u64 delta = a * (2ULL * y_next);
            d += delta;
            d4 += 4ULL * delta;
        }
    }

    return total;
}

u128 psi_prefix_exact(const u64 m) {
    if (m < 3ULL) {
        return 0;
    }
    return easy_sum(m) + 3ULL * hard_one_corner(m);
}

u128 psi_exact(const u64 n) {
    if (n < 3ULL) {
        return 0;
    }
    return psi_prefix_exact(n) - psi_prefix_exact(n - 1ULL);
}

u64 psi_prefix_mod(const u64 m) {
    return static_cast<u64>(psi_prefix_exact(m) % kMod);
}

}  // namespace

int main() {
    assert(psi_exact(3ULL) == 7ULL);
    assert(psi_exact(6ULL) == 34ULL);
    assert(psi_exact(10ULL) == 90ULL);

    assert(psi_prefix_exact(10ULL) == 345ULL);
    assert(psi_prefix_exact(1000ULL) == 172'166'601ULL);

    std::cout << psi_prefix_mod(100'000'000ULL) << '\n';
    return 0;
}
