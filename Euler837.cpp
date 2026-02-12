#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kMod = 1'234'567'891ULL;
static constexpr u64 kInv3 = 823'045'261ULL;
static constexpr u64 kBlock = 2'000'000ULL;

static inline u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

static u64 mod_pow(u64 base, u64 exp) {
    u64 out = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) out = mod_mul(out, base);
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return out;
}

static inline u64 mod_inv(u64 x) {
    return mod_pow(x, kMod - 2);
}

static std::vector<u64> inverse_range(u64 start, u64 len) {
    std::vector<u64> pref(len);
    pref[0] = start % kMod;
    for (u64 i = 1; i < len; ++i) {
        pref[i] = mod_mul(pref[i - 1], (start + i) % kMod);
    }

    u64 inv_prod = mod_inv(pref[len - 1]);
    for (u64 i = len; i-- > 0;) {
        const u64 prev = (i == 0) ? 1 : pref[i - 1];
        const u64 v = (start + i) % kMod;
        pref[i] = mod_mul(inv_prod, prev);
        inv_prod = mod_mul(inv_prod, v);
    }
    return pref;
}

static u64 binom_mod(u64 n, u64 k) {
    if (k > n) return 0;
    k = std::min(k, n - k);
    if (k == 0) return 1;

    u64 out = 1;
    const u64 offset = n - k;
    for (u64 l = 1; l <= k; l += kBlock) {
        const u64 r = std::min(k, l + kBlock - 1);
        const u64 len = r - l + 1;
        std::vector<u64> invs = inverse_range(l, len);
        for (u64 i = 0; i < len; ++i) {
            out = mod_mul(out, (offset + l + i) % kMod);
            out = mod_mul(out, invs[i]);
        }
    }
    return out;
}

static u64 coeff_q(u64 m, u64 n) {
    const u64 a = std::min(m, n);
    const u64 b = std::max(m, n);
    const u64 total = a + b;
    if (total & 1ULL) return 0;
    const u64 t = total >> 1;
    const u64 max_k = a >> 1;

    u64 term = binom_mod(t, a);
    u64 sum = term;

    if (max_k > 0) {
        const u64 base2 = t - a + 1;
        for (u64 l = 0; l < max_k; l += kBlock) {
            const u64 r = std::min(max_k - 1, l + kBlock - 1);
            const u64 len = r - l + 1;
            std::vector<u64> inv1 = inverse_range(l + 1, len);
            std::vector<u64> inv2 = inverse_range(base2 + l, len);
            for (u64 i = 0; i < len; ++i) {
                const u64 k = l + i;
                const u64 u = a - 2 * k;
                term = mod_mul(term, u % kMod);
                term = mod_mul(term, (u - 1) % kMod);
                term = mod_mul(term, inv1[i]);
                term = mod_mul(term, inv2[i]);
                sum += term;
                if (sum >= kMod) sum -= kMod;
            }
        }
    }

    return (a & 1ULL) ? (sum == 0 ? 0 : kMod - sum) : sum;
}

static u64 solve(u64 m, u64 n) {
    if ((m + n) & 1ULL) return 0;
    const u64 c = coeff_q(m, n);
    const u64 b = binom_mod(m + n, m);
    return mod_mul((b + mod_mul(2, c)) % kMod, kInv3);
}

static u64 brute_small(int m, int n) {
    const std::array<std::array<int, 3>, 6> perms = {{
        {{0, 1, 2}}, {{0, 2, 1}}, {{1, 0, 2}},
        {{1, 2, 0}}, {{2, 0, 1}}, {{2, 1, 0}}
    }};

    std::array<int, 6> trans_ab{};
    std::array<int, 6> trans_bc{};
    for (int i = 0; i < 6; ++i) {
        auto p1 = perms[i];
        std::swap(p1[0], p1[1]);
        auto p2 = perms[i];
        std::swap(p2[1], p2[2]);
        for (int j = 0; j < 6; ++j) {
            if (perms[j] == p1) trans_ab[i] = j;
            if (perms[j] == p2) trans_bc[i] = j;
        }
    }

    std::vector<std::vector<std::array<u64, 6>>> dp(
        m + 1, std::vector<std::array<u64, 6>>(n + 1, {0, 0, 0, 0, 0, 0}));
    dp[0][0][0] = 1;

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0 && j == 0) continue;
            std::array<u64, 6> cur = {0, 0, 0, 0, 0, 0};
            if (i > 0) {
                for (int s = 0; s < 6; ++s) {
                    cur[trans_ab[s]] += dp[i - 1][j][s];
                }
            }
            if (j > 0) {
                for (int s = 0; s < 6; ++s) {
                    cur[trans_bc[s]] += dp[i][j - 1][s];
                }
            }
            dp[i][j] = cur;
        }
    }
    return dp[m][n][0];
}

int main() {
    for (int m = 0; m <= 8; ++m) {
        for (int n = 0; n <= 8; ++n) {
            assert(solve(m, n) == brute_small(m, n) % kMod);
        }
    }
    assert(solve(3, 3) == 2);
    assert(solve(123, 321) == 172'633'303ULL);

    std::cout << solve(123'456'789ULL, 987'654'321ULL) << '\n';
    return 0;
}
