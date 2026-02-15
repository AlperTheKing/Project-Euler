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
static constexpr u64 kBlock = 8'000'000ULL;
static constexpr u64 kBarrett = static_cast<u64>((static_cast<u128>(1) << 64) / kMod);

static inline u64 mod_reduce(u64 x) {
    const u64 q = static_cast<u64>((static_cast<u128>(x) * kBarrett) >> 64);
    u64 r = x - q * kMod;
    if (r >= kMod) r -= kMod;
    if (r >= kMod) r -= kMod;
    return r;
}

static inline u64 mod_mul(u64 a, u64 b) {
    return mod_reduce(a * b);
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

static void inverse_range(u64 start, u64 len, std::vector<u64>& pref) {
    pref.resize(static_cast<std::size_t>(len));
    u64 cur = start % kMod;
    pref[0] = cur;
    for (u64 i = 1; i < len; ++i) {
        if (++cur == kMod) cur = 0;
        pref[static_cast<std::size_t>(i)] = mod_mul(pref[static_cast<std::size_t>(i - 1)], cur);
    }

    u64 inv_prod = mod_inv(pref[static_cast<std::size_t>(len - 1)]);
    cur = (start + len - 1) % kMod;
    for (u64 i = len; i-- > 0;) {
        const u64 prev = (i == 0) ? 1 : pref[static_cast<std::size_t>(i - 1)];
        pref[static_cast<std::size_t>(i)] = mod_mul(inv_prod, prev);
        inv_prod = mod_mul(inv_prod, cur);
        if (cur == 0) cur = kMod - 1;
        else --cur;
    }
}

static u64 binom_mod(u64 n, u64 k) {
    if (k > n) return 0;
    k = std::min(k, n - k);
    if (k == 0) return 1;

    u64 out = 1;
    const u64 offset = n - k;
    std::vector<u64> invs;
    invs.reserve(static_cast<std::size_t>(kBlock));
    for (u64 l = 1; l <= k; l += kBlock) {
        const u64 r = std::min(k, l + kBlock - 1);
        const u64 len = r - l + 1;
        inverse_range(l, len, invs);
        u64 num = (offset + l) % kMod;
        for (u64 i = 0; i < len; ++i) {
            out = mod_mul(out, num);
            out = mod_mul(out, invs[static_cast<std::size_t>(i)]);
            if (++num == kMod) num = 0;
        }
    }
    return out;
}

static u64 solve(u64 m, u64 n) {
    if ((m + n) & 1ULL) return 0;

    const u64 t = (m + n) >> 1;
    const u64 parity = m & 1ULL;
    const u64 max_k = std::min(m, n);

    u64 k = parity;
    u64 a = (m - parity) >> 1;
    u64 b = (n - parity) >> 1;

    u64 weight = 0;
    if (parity == 0ULL) {
        weight = binom_mod(t, a);
    } else {
        weight = mod_mul(t % kMod, binom_mod(t - 1ULL, a));
    }

    u64 pow2k = (parity == 0ULL) ? 1ULL : 2ULL;
    const u64 sign_part = (parity == 0ULL) ? 2ULL : (kMod - 2ULL);  // 2*(-1)^k

    u64 answer = 0ULL;
    auto add_current = [&]() {
        u64 f3 = pow2k + sign_part;
        if (f3 >= kMod) f3 -= kMod;
        f3 = mod_mul(f3, kInv3);
        const u64 inc = mod_mul(weight, f3);
        answer += inc;
        if (answer >= kMod) answer -= kMod;
    };

    add_current();
    if (k == max_k) return answer;

    u64 updates = (max_k - k) >> 1;
    std::vector<u64> invs;
    invs.reserve(static_cast<std::size_t>(2 * kBlock));
    u64 a_mod = a % kMod;
    u64 b_mod = b % kMod;

    while (updates > 0ULL) {
        const u64 len = std::min<u64>(updates, kBlock);
        inverse_range(k + 1ULL, 2ULL * len, invs);

        for (u64 i = 0; i < len; ++i) {
            const u64 inv_k1 = invs[static_cast<std::size_t>(2ULL * i)];
            const u64 inv_k2 = invs[static_cast<std::size_t>(2ULL * i + 1ULL)];

            weight = mod_mul(weight, a_mod);
            weight = mod_mul(weight, b_mod);
            weight = mod_mul(weight, inv_k1);
            weight = mod_mul(weight, inv_k2);

            --a;
            --b;
            if (a_mod == 0) a_mod = kMod - 1;
            else --a_mod;
            if (b_mod == 0) b_mod = kMod - 1;
            else --b_mod;
            k += 2ULL;
            pow2k = mod_mul(pow2k, 4ULL);

            add_current();
        }
        updates -= len;
    }

    return answer;
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
