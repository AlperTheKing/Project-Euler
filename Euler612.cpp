#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 612: count pairs of integers < 10^K whose decimal digit-sets intersect.

using i64 = long long;
using i128 = __int128_t;

static constexpr i64 MOD = 1000267129LL;
static constexpr int DIGS = 10;
static constexpr int FULL = (1 << DIGS) - 1;

static inline i64 mod_norm(i64 x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return x;
}

static inline i64 mul_mod(i64 a, i64 b) { return (i64)((i128)a * b % MOD); }

static i64 pow_mod(i64 a, std::uint64_t e) {
    i64 r = 1 % MOD;
    a %= MOD;
    while (e) {
        if (e & 1) r = mul_mod(r, a);
        a = mul_mod(a, a);
        e >>= 1;
    }
    return r;
}

static std::vector<i64> count_by_mask(int max_len) {
    std::vector<i64> cnt(1 << DIGS, 0);
    std::vector<i64> dp(1 << DIGS, 0), ndp(1 << DIGS, 0);

    for (int d = 1; d <= 9; ++d) {
        dp[1 << d] += 1;
        if (dp[1 << d] >= MOD) dp[1 << d] -= MOD;
    }

    for (int len = 1; len <= max_len; ++len) {
        for (int m = 0; m <= FULL; ++m) {
            cnt[m] += dp[m];
            if (cnt[m] >= MOD) cnt[m] -= MOD;
        }
        if (len == max_len) break;

        std::fill(ndp.begin(), ndp.end(), 0);
        for (int m = 0; m <= FULL; ++m) {
            const i64 v = dp[m];
            if (!v) continue;
            for (int d = 0; d <= 9; ++d) {
                const int nm = m | (1 << d);
                ndp[nm] += v;
                if (ndp[nm] >= MOD) ndp[nm] -= MOD;
            }
        }
        dp.swap(ndp);
    }
    return cnt;
}

static i64 f_power10(int K) {
    const i64 inv2 = (MOD + 1) / 2;

    const std::vector<i64> cnt = count_by_mask(K);

    std::vector<i64> sum_sub = cnt; // SOS DP: sum_sub[mask] = sum_{s subset mask} cnt[s]
    for (int b = 0; b < DIGS; ++b) {
        for (int mask = 0; mask <= FULL; ++mask) {
            if (mask & (1 << b)) {
                sum_sub[mask] += sum_sub[mask ^ (1 << b)];
                if (sum_sub[mask] >= MOD) sum_sub[mask] -= MOD;
            }
        }
    }

    i64 ordered_disjoint = 0; // ordered pairs (p,q) with disjoint digit-sets
    for (int a = 0; a <= FULL; ++a) {
        ordered_disjoint += mul_mod(cnt[a], sum_sub[FULL ^ a]);
        ordered_disjoint %= MOD;
    }
    const i64 disjoint_pairs = mul_mod(ordered_disjoint, inv2);

    const i64 M = mod_norm(pow_mod(10, (std::uint64_t)K) - 1); // count of integers in [1,10^K-1]
    const i64 total_pairs = mul_mod(mul_mod(M, mod_norm(M - 1)), inv2);

    return mod_norm(total_pairs - disjoint_pairs);
}

static int digit_mask_u64(std::uint64_t x) {
    int m = 0;
    while (x) {
        m |= 1 << (int)(x % 10);
        x /= 10;
    }
    return m;
}

static std::uint64_t brute_f(int n) {
    std::vector<int> masks((std::size_t)n, 0);
    for (int i = 1; i < n; ++i) masks[i] = digit_mask_u64((std::uint64_t)i);

    std::uint64_t ans = 0;
    for (int p = 1; p < n; ++p) {
        for (int q = p + 1; q < n; ++q) {
            if (masks[p] & masks[q]) ++ans;
        }
    }
    return ans;
}

int main() {
    assert(f_power10(2) == 1539);
    assert(f_power10(3) == (i64)(brute_f(1000) % (std::uint64_t)MOD));

    std::cout << f_power10(18) << "\n";
    return 0;
}

