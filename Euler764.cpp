#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1'000'000'000ULL;

u128 pow4(u64 x) {
    const u128 a = static_cast<u128>(x);
    return a * a * a * a;
}

u64 iroot4(u128 n) {
    if (n == 0) {
        return 0;
    }
    long double d = static_cast<long double>(n);
    u64 x = static_cast<u64>(std::pow(d, 0.25L));
    while (pow4(x + 1) <= n) {
        ++x;
    }
    while (pow4(x) > n) {
        --x;
    }
    return x;
}

u64 isqrt_u128(u128 n) {
    if (n == 0) {
        return 0;
    }
    long double d = static_cast<long double>(n);
    u64 x = static_cast<u64>(std::sqrt(d));
    while (static_cast<u128>(x + 1) * (x + 1) <= n) {
        ++x;
    }
    while (static_cast<u128>(x) * x > n) {
        --x;
    }
    return x;
}

template <typename F>
void enumerate_solutions(u64 n, F&& add_term) {
    const u64 q_max = iroot4(static_cast<u128>(2) * n);
    for (u64 q = 1; q <= q_max; q += 2) {
        const u128 q4 = pow4(q);
        if (q4 > static_cast<u128>(2) * n) {
            continue;
        }
        const u128 lim = static_cast<u128>(2) * n - q4;
        u64 p_max = iroot4(lim);
        if (p_max >= q) {
            p_max = q - 1;
        }
        if ((p_max & 1ULL) == 0ULL) {
            if (p_max == 0) {
                continue;
            }
            --p_max;
        }

        for (u64 p = 1; p <= p_max; p += 2) {
            if (std::gcd(p, q) != 1ULL) {
                continue;
            }
            const u128 p4 = pow4(p);
            const u128 z = (p4 + q4) / 2;
            if (z > n) {
                continue;
            }
            const u128 x = (q4 - p4) / 8;
            const u128 y = static_cast<u128>(p) * q;
            if (y > n) {
                continue;
            }
            add_term(x + y + z);
        }
    }

    const u64 s_max = iroot4(n / 4);
    for (u64 s = 1; s <= s_max; s += 2) {
        const u128 s4 = pow4(s);
        if (s4 > static_cast<u128>(n / 4)) {
            break;
        }
        const u128 rem = static_cast<u128>(n / 4) - s4;
        const u64 r_max = iroot4(rem / 4);

        for (u64 r = 1; r <= r_max; ++r) {
            if (std::gcd(r, s) != 1ULL) {
                continue;
            }
            const u128 r4 = pow4(r);
            const u128 z = 4 * (s4 + 4 * r4);
            if (z > n) {
                continue;
            }
            const u128 a = 4 * r4;
            const u128 x = (s4 >= a) ? (s4 - a) : (a - s4);
            const u128 y = static_cast<u128>(4) * r * s;
            if (y > n) {
                continue;
            }
            add_term(x + y + z);
        }
    }
}

u128 solve_exact(u64 n) {
    u128 total = 0;
    enumerate_solutions(n, [&](u128 term) { total += term; });
    return total;
}

u64 solve_mod(u64 n, u64 mod) {
    u64 total = 0;
    enumerate_solutions(n, [&](u128 term) {
        total += static_cast<u64>(term % mod);
        if (total >= mod) {
            total %= mod;
        }
    });
    return total % mod;
}

u64 brute_small(u64 n) {
    u64 total = 0;
    const u64 y_max = static_cast<u64>(std::sqrt(static_cast<long double>(n)));

    for (u64 y = 1; y <= y_max; ++y) {
        const u128 y2 = static_cast<u128>(y) * y;
        const u128 y4 = y2 * y2;
        for (u64 x = 1; x <= n; ++x) {
            const u128 z2 = static_cast<u128>(16) * x * x + y4;
            const u64 z = isqrt_u128(z2);
            if (z > n) {
                continue;
            }
            if (static_cast<u128>(z) * z != z2) {
                continue;
            }
            if (std::gcd(std::gcd(x, y), z) != 1ULL) {
                continue;
            }
            total += x + y + z;
        }
    }
    return total;
}

}  // namespace

int main() {
    assert(static_cast<u64>(solve_exact(100)) == 81ULL);
    assert(static_cast<u64>(solve_exact(10'000)) == 112'851ULL);
    assert(solve_mod(10'000'000, MOD) == 248'876'211ULL);
    assert(static_cast<u64>(solve_exact(200)) == brute_small(200));

    std::cout << solve_mod(10'000'000'000'000'000ULL, MOD) << '\n';
    return 0;
}
