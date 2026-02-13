#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

u64 isqrt_u128(u128 x) {
    long double y = std::sqrt(static_cast<long double>(x));
    u64 r = static_cast<u64>(y);
    while ((u128)(r + 1) * (r + 1) <= x) {
        ++r;
    }
    while ((u128)r * r > x) {
        --r;
    }
    return r;
}

u64 best_for_n(u64 R, u64 n) {
    const u128 C = (u128)2 * R;
    const u128 nn = (u128)n * n;
    if (nn >= C) {
        return 0;
    }

    const u128 t = C - 1 - nn;
    u64 m = isqrt_u128(t);
    if (m <= n) {
        return 0;
    }

    if (((m ^ n) & 1ULL) == 0ULL) {
        --m;
    }

    while (m > n && std::gcd(m, n) != 1) {
        m -= 2;
    }

    if (m <= n) {
        return 0;
    }

    return n * (m - n);
}

u64 exhaustive_F(u64 R) {
    const u128 C = (u128)2 * R;
    const u64 n_max = isqrt_u128((C - 1) / 2);

    u64 best = 0;
    for (u64 n = 1; n <= n_max; ++n) {
        const u64 cand = best_for_n(R, n);
        if (cand > best) {
            best = cand;
        }
    }
    return best;
}

u64 fast_F(u64 R) {
    const u128 C = (u128)2 * R;
    const long double sqrt2 = std::sqrt(2.0L);
    const long double denom = 4.0L + 2.0L * sqrt2;
    const long double Cld = static_cast<long double>(C - 1);

    u64 n0 = static_cast<u64>(std::llround(std::sqrt(Cld / denom)));
    if (n0 == 0) {
        n0 = 1;
    }

    u64 best = 0;
    constexpr int kWindow = 4096;
    for (int d = -kWindow; d <= kWindow; ++d) {
        const std::int64_t ni = static_cast<std::int64_t>(n0) + d;
        if (ni <= 0) {
            continue;
        }
        const u64 cand = best_for_n(R, static_cast<u64>(ni));
        if (cand > best) {
            best = cand;
        }
    }

    const auto ub_real = [&](long double n) {
        if (n <= 0) {
            return 0.0L;
        }
        const long double t = Cld - n * n;
        if (t <= 0) {
            return 0.0L;
        }
        return n * (std::sqrt(t) - n);
    };

    const long double n0ld = static_cast<long double>(n0);
    const long double nmaxld = std::sqrt(Cld / 2.0L);

    long double left = n0ld;
    long double right = n0ld;

    if (ub_real(n0ld) > static_cast<long double>(best)) {
        long double lo = 0.0L;
        long double hi = n0ld;
        for (int it = 0; it < 140; ++it) {
            const long double mid = (lo + hi) * 0.5L;
            if (ub_real(mid) > static_cast<long double>(best)) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        left = hi;

        lo = n0ld;
        hi = nmaxld;
        for (int it = 0; it < 140; ++it) {
            const long double mid = (lo + hi) * 0.5L;
            if (ub_real(mid) > static_cast<long double>(best)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        right = lo;
    }

    u64 L = 1;
    if (left > 16.0L) {
        L = static_cast<u64>(std::floor(left)) - 16;
    }

    const u64 n_max = isqrt_u128((C - 1) / 2);
    u64 Rn = static_cast<u64>(std::ceil(right)) + 16;
    if (Rn > n_max) {
        Rn = n_max;
    }

    for (u64 n = L; n <= Rn; ++n) {
        const u128 nn = (u128)n * n;
        if (nn >= C) {
            break;
        }
        const u64 mmax = isqrt_u128(C - 1 - nn);
        if (mmax <= n) {
            continue;
        }
        const u64 ub = n * (mmax - n);
        if (ub <= best) {
            continue;
        }

        const u64 cand = best_for_n(R, n);
        if (cand > best) {
            best = cand;
        }
    }

    return best;
}

void validate() {
    assert(exhaustive_F(100) == 36);
    assert(fast_F(100) == 36);

    for (u64 R = 2; R <= 5000; ++R) {
        assert(fast_F(R) == exhaustive_F(R));
    }

    u64 seed = 1;
    for (int i = 0; i < 200; ++i) {
        seed = seed * 6364136223846793005ULL + 1ULL;
        const u64 R = 2 + (seed % 100'000'000ULL);
        assert(fast_F(R) == exhaustive_F(R));
    }
}

}  // namespace

int main() {
    validate();
    std::cout << fast_F(1'000'000'000'000'000'000ULL) << '\n';
    return 0;
}
