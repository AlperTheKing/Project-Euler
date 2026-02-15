#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'000ULL;

u64 isqrt_u64(u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) {
        ++r;
    }
    while (r * r > n) {
        --r;
    }
    return r;
}

u128 triangular_minus_one(u64 n) {
    return static_cast<u128>(n) * static_cast<u128>(n + 1ULL) / 2U - 1U;
}

u128 min_factor_sum(u64 n) {
    const u64 v = isqrt_u64(n);

    std::vector<u64> s_cnt(v + 1ULL, 0ULL);
    std::vector<u128> s_sum(v + 1ULL, 0ULL);
    std::vector<u64> l_cnt(v + 1ULL, 0ULL);
    std::vector<u128> l_sum(v + 1ULL, 0ULL);
    std::vector<std::uint8_t> used(v + 1ULL, 0U);

    for (u64 i = 0; i <= v; ++i) {
        s_cnt[static_cast<std::size_t>(i)] = (i == 0ULL) ? 0ULL : (i - 1ULL);
        s_sum[static_cast<std::size_t>(i)] = triangular_minus_one(i);
        if (i == 0ULL) {
            l_cnt[0] = 0ULL;
            l_sum[0] = 0U;
        } else {
            const u64 q = n / i;
            l_cnt[static_cast<std::size_t>(i)] = q - 1ULL;
            l_sum[static_cast<std::size_t>(i)] = triangular_minus_one(q);
        }
    }

    u128 ret = 0U;
    for (u64 p = 2ULL; p <= v; ++p) {
        const std::size_t ps = static_cast<std::size_t>(p);
        if (s_cnt[ps] == s_cnt[ps - 1ULL]) {
            continue;
        }

        const u64 p_cnt = s_cnt[ps - 1ULL];
        const u128 p_sum = s_sum[ps - 1ULL];
        const u64 q = p * p;

        ret += static_cast<u128>(p) * static_cast<u128>(l_cnt[ps] - p_cnt);

        l_cnt[1] -= (l_cnt[ps] - p_cnt);
        l_sum[1] -= (l_sum[ps] - p_sum) * static_cast<u128>(p);

        const u64 interval = (p & 1ULL) + 1ULL;
        const u64 end = std::min(v, n / q);
        for (u64 i = p + interval; i <= end; i += interval) {
            const std::size_t is = static_cast<std::size_t>(i);
            if (used[is] != 0U) {
                continue;
            }
            const u64 d = i * p;
            if (d <= v) {
                const std::size_t ds = static_cast<std::size_t>(d);
                l_cnt[is] -= (l_cnt[ds] - p_cnt);
                l_sum[is] -= (l_sum[ds] - p_sum) * static_cast<u128>(p);
            } else {
                const u64 t = n / d;
                const std::size_t ts = static_cast<std::size_t>(t);
                l_cnt[is] -= (s_cnt[ts] - p_cnt);
                l_sum[is] -= (s_sum[ts] - p_sum) * static_cast<u128>(p);
            }
        }

        if (q <= v) {
            const u64 step = p * interval;
            for (u64 i = q; i < end; i += step) {
                used[static_cast<std::size_t>(i)] = 1U;
            }
        }

        for (u64 i = v; i >= q; --i) {
            const std::size_t is = static_cast<std::size_t>(i);
            const std::size_t ts = static_cast<std::size_t>(i / p);
            s_cnt[is] -= (s_cnt[ts] - p_cnt);
            s_sum[is] -= (s_sum[ts] - p_sum) * static_cast<u128>(p);
        }
    }

    return l_sum[1] + ret;
}

u64 brute_small(int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            if (static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(n)) {
                for (int m = i * i; m <= n; m += i) {
                    if (spf[static_cast<std::size_t>(m)] == 0) {
                        spf[static_cast<std::size_t>(m)] = i;
                    }
                }
            }
        }
    }

    u64 sum = 0ULL;
    for (int i = 2; i <= n; ++i) {
        sum += static_cast<u64>(spf[static_cast<std::size_t>(i)]);
    }
    return sum;
}

u64 mod_u128(u128 x, u64 mod) {
    return static_cast<u64>(x % static_cast<u128>(mod));
}

bool run_checkpoints() {
    if (min_factor_sum(100ULL) != 1'257ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    if (min_factor_sum(100ULL) != static_cast<u128>(brute_small(100))) {
        std::cerr << "Checkpoint failed: solve/brute mismatch at 100\n";
        return false;
    }
    if (min_factor_sum(1'000ULL) != static_cast<u128>(brute_small(1'000))) {
        std::cerr << "Checkpoint failed: solve/brute mismatch at 1000\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 1'000'000'000'000ULL;
    const u128 full = min_factor_sum(n);
    const u64 answer = mod_u128(full, kMod);
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}
