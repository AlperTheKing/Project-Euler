#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

int threshold(const int n) {
    return (n + 1) * (n + 2) / 2;
}

int max_n_for_m(const int m) {
    int n = 0;
    while (threshold(n + 1) <= m) {
        ++n;
    }
    return n;
}

template <bool UseMod>
u64 solve_a2(const int target_m, const u64 mod_value = 0ULL) {
    const int nmax = max_n_for_m(target_m);
    const int ring = nmax + 8;
    const int nalloc = nmax + 2;

    std::vector<std::vector<u64>> u(nalloc + 1);
    std::vector<std::vector<u64>> v(nalloc + 1);
    for (int n = 0; n <= nalloc; ++n) {
        u[n].assign(static_cast<std::size_t>(n + 2) * static_cast<std::size_t>(ring), 0ULL);
        v[n].assign(static_cast<std::size_t>(n + 2) * static_cast<std::size_t>(ring), 0ULL);
    }

    std::vector<u64> a(static_cast<std::size_t>(ring), 0ULL);
    std::vector<u64> f0(static_cast<std::size_t>(ring), 0ULL);
    a[0] = 1ULL;

    auto norm = [&](const u128 x) -> u64 {
        if constexpr (UseMod) {
            return static_cast<u64>(x % static_cast<u128>(mod_value));
        } else {
            return static_cast<u64>(x);
        }
    };

    auto get_a = [&](const int m) -> u64 {
        if (m < 0) {
            return 0ULL;
        }
        return a[static_cast<std::size_t>(m % ring)];
    };

    auto get_f0 = [&](const int m) -> u64 {
        if (m < 0) {
            return 0ULL;
        }
        return f0[static_cast<std::size_t>(m % ring)];
    };

    auto get_u = [&](const int n, const int k, const int m) -> u64 {
        if (m < 0) {
            return 0ULL;
        }
        if (n == 0) {
            return (k == 0) ? get_f0(m) : 0ULL;
        }
        if (n < 0 || n > nalloc || k < 1 || k > n || m < threshold(n)) {
            return 0ULL;
        }
        const std::size_t idx = static_cast<std::size_t>(k) * static_cast<std::size_t>(ring)
                              + static_cast<std::size_t>(m % ring);
        return u[static_cast<std::size_t>(n)][idx];
    };

    auto get_v = [&](const int n, const int k, const int m) -> u64 {
        if (m < 0) {
            return 0ULL;
        }
        if (n == 0) {
            return (k == 0) ? get_f0(m) : 0ULL;
        }
        if (n < 0 || n > nalloc || k < 1 || k > n || m < threshold(n)) {
            return 0ULL;
        }
        const std::size_t idx = static_cast<std::size_t>(k) * static_cast<std::size_t>(ring)
                              + static_cast<std::size_t>(m % ring);
        return v[static_cast<std::size_t>(n)][idx];
    };

    for (int m = 0; m <= target_m; ++m) {
        const int cur = m % ring;

        for (int n = 1; n <= nmax; ++n) {
            for (int k = 1; k <= n; ++k) {
                const std::size_t idx = static_cast<std::size_t>(k) * static_cast<std::size_t>(ring)
                                      + static_cast<std::size_t>(cur);
                u[static_cast<std::size_t>(n)][idx] = 0ULL;
                v[static_cast<std::size_t>(n)][idx] = 0ULL;
            }
        }

        for (int n = 1; n <= nmax; ++n) {
            if (m < threshold(n)) {
                continue;
            }

            for (int k = 1; k <= n; ++k) {
                const int h = (k < n) ? k : (k - 1);

                const u128 u_sum = static_cast<u128>(get_u(n, k, m - n - 2))
                                 + static_cast<u128>(get_v(n + 1, 1, m - n - 3))
                                 + static_cast<u128>(get_u(n + 1, k + 1, m - n - 3))
                                 + static_cast<u128>(get_u(n - 1, h, m - n - 1))
                                 + static_cast<u128>(get_v(n, 1, m - n - 2))
                                 + static_cast<u128>(get_u(n, h + 1, m - n - 2));

                const u128 v_sum = static_cast<u128>(get_v(n, k, m - n - 2))
                                 + static_cast<u128>(get_v(n + 1, k + 1, m - n - 3))
                                 + static_cast<u128>((k == n) ? 2ULL : 1ULL)
                                       * static_cast<u128>(get_u(n + 1, 1, m - n - 3))
                                 + static_cast<u128>(get_v(n - 1, h, m - n - 1))
                                 + static_cast<u128>(get_v(n, h + 1, m - n - 2))
                                 + static_cast<u128>((k == n) ? 2ULL : 1ULL)
                                       * static_cast<u128>(get_u(n, 1, m - n - 2));

                const std::size_t idx = static_cast<std::size_t>(k) * static_cast<std::size_t>(ring)
                                      + static_cast<std::size_t>(cur);
                u[static_cast<std::size_t>(n)][idx] = norm(u_sum);
                v[static_cast<std::size_t>(n)][idx] = norm(v_sum);
            }
        }

        const u128 f0_sum = static_cast<u128>(get_a(m - 1))
                          + static_cast<u128>(4ULL) * static_cast<u128>(get_f0(m - 2))
                          + static_cast<u128>(2ULL) * static_cast<u128>(get_u(1, 1, m - 3))
                          + static_cast<u128>(get_v(1, 1, m - 3));
        f0[static_cast<std::size_t>(cur)] = norm(f0_sum);

        if (m >= 1) {
            const u128 a_sum = static_cast<u128>(3ULL) * static_cast<u128>(get_a(m - 1))
                             + static_cast<u128>(3ULL) * static_cast<u128>(get_f0(m - 2));
            a[static_cast<std::size_t>(cur)] = norm(a_sum);
        }
    }

    return a[static_cast<std::size_t>(target_m % ring)];
}

u64 D_exact(const int n) {
    assert(n >= 1);
    return solve_a2<false>(n - 1);
}

u64 D_mod(const int n, const u64 mod) {
    assert(n >= 1);
    return solve_a2<true>(n - 1, mod);
}

}  // namespace

int main() {
    constexpr u64 MOD = 1'000'000'000ULL;

    assert(D_exact(10) == 44'499ULL);
    assert(D_exact(20) == 9'204'559'704ULL);
    assert(D_exact(32) == 22'037'102'049'132'222ULL);
    assert(D_mod(100, MOD) == 780'166'455ULL);

    const u64 ans = D_mod(10'000, MOD);
    std::cout << std::setw(9) << std::setfill('0') << ans << '\n';
    return 0;
}
