#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr i64 kMod = 1'000'000'000LL;

inline i64 mod_norm(i64 x) {
    x %= kMod;
    if (x < 0) x += kMod;
    return x;
}

inline i64 mod_add(i64 a, i64 b) {
    a += b;
    if (a >= kMod) a -= kMod;
    return a;
}

inline i64 mod_sub(i64 a, i64 b) {
    a -= b;
    if (a < 0) a += kMod;
    return a;
}

inline i64 mod_mul(i64 a, i64 b) {
    return static_cast<i64>((static_cast<__int128>(a) * b) % kMod);
}

u64 isqrt_u64(u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

inline i64 sum_2_to_mod(u64 x) {
    if (x < 2ULL) return 0;
    const u128 s = static_cast<u128>(x) * static_cast<u128>(x + 1ULL) / 2U;
    return static_cast<i64>((s - 1U) % static_cast<u128>(kMod));
}

struct Solver {
    u64 n = 0;
    u64 m = 0;

    std::vector<i64> pre_cnt;
    std::vector<i64> pre_sum;
    std::vector<i64> hou_cnt;
    std::vector<i64> hou_sum;
    std::vector<int> primes;

    i64 dfs(u64 res, int last, int from_prime) {
        i64 ret = 0;
        if (from_prime > 0) {
            const i64 base = (res > m) ? hou_sum[static_cast<std::size_t>(n / res)]
                                       : pre_sum[static_cast<std::size_t>(res)];
            ret = mod_sub(base, pre_sum[static_cast<std::size_t>(primes[static_cast<std::size_t>(last)] - 1)]);
        } else {
            ret = hou_sum[1];
        }

        for (int i = last; i < static_cast<int>(primes.size()); ++i) {
            const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(i)]);
            if (p * p > res) break;

            u64 nres = res;
            for (u64 q = p; q * p <= res; q *= p) {
                nres /= p;
                ret = mod_add(ret, dfs(nres, i + 1, static_cast<int>(p)));
                ret = mod_add(ret, static_cast<i64>(p % static_cast<u64>(kMod)));
            }
        }
        return ret;
    }

    i64 solve(u64 input_n) {
        n = input_n;
        if (n <= 1ULL) return 0;
        m = isqrt_u64(n);

        pre_cnt.assign(static_cast<std::size_t>(m + 1), 0);
        pre_sum.assign(static_cast<std::size_t>(m + 1), 0);
        hou_cnt.assign(static_cast<std::size_t>(m + 1), 0);
        hou_sum.assign(static_cast<std::size_t>(m + 1), 0);
        primes.clear();
        primes.reserve(static_cast<std::size_t>(m / 10 + 16));

        for (u64 i = 1; i <= m; ++i) {
            pre_cnt[static_cast<std::size_t>(i)] = static_cast<i64>(i - 1ULL);
            pre_sum[static_cast<std::size_t>(i)] = sum_2_to_mod(i);
            const u64 v = n / i;
            hou_cnt[static_cast<std::size_t>(i)] = static_cast<i64>(v - 1ULL);
            hou_sum[static_cast<std::size_t>(i)] = sum_2_to_mod(v);
        }

        for (u64 p = 2; p <= m; ++p) {
            const std::size_t ps = static_cast<std::size_t>(p);
            if (pre_cnt[ps] == pre_cnt[ps - 1]) continue;

            primes.push_back(static_cast<int>(p));

            const u64 p2 = p * p;
            const u64 q = n / p;
            const i64 p_cnt = pre_cnt[ps - 1];
            const i64 p_sum = pre_sum[ps - 1];
            const u64 mid = m / p;
            const u64 end = std::min(m, n / p2);
            const i64 p_mod = static_cast<i64>(p % static_cast<u64>(kMod));

            for (u64 i = 1; i <= mid; ++i) {
                const std::size_t is = static_cast<std::size_t>(i);
                const std::size_t ips = static_cast<std::size_t>(i * p);
                hou_cnt[is] -= (hou_cnt[ips] - p_cnt);
                hou_sum[is] = mod_sub(hou_sum[is], mod_mul(mod_sub(hou_sum[ips], p_sum), p_mod));
            }
            for (u64 i = mid + 1; i <= end; ++i) {
                const std::size_t is = static_cast<std::size_t>(i);
                const std::size_t qi = static_cast<std::size_t>(q / i);
                hou_cnt[is] -= (pre_cnt[qi] - p_cnt);
                hou_sum[is] = mod_sub(hou_sum[is], mod_mul(mod_sub(pre_sum[qi], p_sum), p_mod));
            }
            for (u64 i = m; i >= p2; --i) {
                const std::size_t is = static_cast<std::size_t>(i);
                const std::size_t ip = static_cast<std::size_t>(i / p);
                pre_cnt[is] -= (pre_cnt[ip] - p_cnt);
                pre_sum[is] = mod_sub(pre_sum[is], mod_mul(mod_sub(pre_sum[ip], p_sum), p_mod));
            }
        }

        primes.push_back(static_cast<int>(m + 1ULL));
        return mod_norm(dfs(n, 0, 0));
    }
};

i64 brute(u64 n) {
    std::vector<u64> lp(static_cast<std::size_t>(n + 1), 0ULL);
    for (u64 i = 2; i <= n; ++i) {
        if (lp[static_cast<std::size_t>(i)] != 0ULL) continue;
        for (u64 j = i; j <= n; j += i) {
            lp[static_cast<std::size_t>(j)] = i;
        }
    }
    i64 s = 0;
    for (u64 i = 2; i <= n; ++i) {
        s = mod_add(s, static_cast<i64>(lp[static_cast<std::size_t>(i)] % static_cast<u64>(kMod)));
    }
    return s;
}

}  // namespace

int main() {
    Solver solver;
    assert(solver.solve(10) == 32);
    assert(solver.solve(100) == 1915);
    assert(solver.solve(10'000) == 10'118'280);
    assert(solver.solve(200'000) == brute(200'000));

    constexpr u64 n = 201'820'182'018ULL;
    std::cout << solver.solve(n) << '\n';
    return 0;
}

