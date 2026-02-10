#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 614: special partitions are distinct parts with no part ≡ 2 (mod 4).

using i64 = long long;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using u32 = std::uint32_t;

static constexpr int MOD = 1'000'000'007;

struct NTTMod {
    u32 mod;
    u32 primitive_root;
};

static constexpr NTTMod P1{998244353u, 3u};     // 2^23 * 119 + 1
static constexpr NTTMod P2{1224736769u, 3u};    // 2^24 * 73 + 1
static constexpr NTTMod P3{469762049u, 3u};     // 2^26 * 7 + 1

static u32 mod_pow_u32(u32 a, u32 e, u32 mod) {
    u64 r = 1, x = a;
    while (e) {
        if (e & 1u) r = (r * x) % mod;
        x = (x * x) % mod;
        e >>= 1u;
    }
    return (u32)r;
}

static u32 mod_inv_u32(u32 a, u32 mod) { return mod_pow_u32(a, mod - 2, mod); }

static void ntt_many(std::vector<u32*>& arrs, int n, bool invert, const NTTMod& p) {
    const u32 mod = p.mod;

    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            for (u32* a : arrs) std::swap(a[i], a[j]);
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        u32 wlen = mod_pow_u32(p.primitive_root, (mod - 1) / (u32)len, mod);
        if (invert) wlen = mod_inv_u32(wlen, mod);

        for (int i = 0; i < n; i += len) {
            u64 w = 1;
            const int half = len >> 1;
            for (int j = 0; j < half; ++j) {
                for (u32* a : arrs) {
                    const u32 u = a[i + j];
                    const u32 v = (u32)(w * a[i + j + half] % mod);
                    u32 x = u + v;
                    if (x >= mod) x -= mod;
                    a[i + j] = x;
                    a[i + j + half] = (u >= v) ? (u - v) : (u + mod - v);
                }
                w = (w * wlen) % mod;
            }
        }
    }

    if (invert) {
        const u32 inv_n = mod_inv_u32((u32)n, mod);
        for (u32* a : arrs) {
            for (int i = 0; i < n; ++i) a[i] = (u32)((u64)a[i] * inv_n % mod);
        }
    }
}

static int crt3_to_mod(u32 a1, u32 a2, u32 a3) {
    static const u64 m1 = P1.mod;
    static const u64 m2 = P2.mod;
    static const u64 m3 = P3.mod;
    static const u64 m1m2 = m1 * m2;
    static const u64 inv_m1_m2 = mod_inv_u32((u32)(m1 % m2), (u32)m2);
    static const u64 inv_m1m2_m3 = mod_inv_u32((u32)(m1m2 % m3), (u32)m3);

    const u64 t1 = (u64)((a2 + m2 - (a1 % m2)) % m2) * inv_m1_m2 % m2;
    const u128 x12 = (u128)a1 + (u128)t1 * m1;  // < m1*m2
    const u64 x12_mod_m3 = (u64)(x12 % m3);
    const u64 t2 = (u64)((a3 + m3 - x12_mod_m3) % m3) * inv_m1m2_m3 % m3;
    const u128 x123 = x12 + (u128)t2 * (u128)m1m2;
    return (int)(x123 % (u64)MOD);
}

struct RingMul {
    std::vector<u32> buf[8];

    static int next_pow2(int x) {
        int n = 1;
        while (n < x) n <<= 1;
        return n;
    }

    std::vector<int> mul_trunc(const std::vector<int>& a, const std::vector<int>& b, int n_trunc) {
        if (n_trunc <= 0) return {};
        const int na = std::min((int)a.size(), n_trunc);
        const int nb = std::min((int)b.size(), n_trunc);
        if (na == 0 || nb == 0) return std::vector<int>(n_trunc, 0);

        if ((i64)na * (i64)nb <= 20'000) {
            std::vector<int> out(n_trunc, 0);
            for (int i = 0; i < na; ++i) {
                if (a[i] == 0) continue;
                const i64 ai = a[i];
                const int maxj = std::min(nb - 1, n_trunc - 1 - i);
                for (int j = 0; j <= maxj; ++j) {
                    if (b[j] == 0) continue;
                    out[i + j] = (out[i + j] + (int)(ai * b[j] % MOD)) % MOD;
                }
            }
            return out;
        }

        const int lenT = (n_trunc + 3) / 4;
        const int L = next_pow2(2 * lenT);

        std::vector<int> A[4], B[4];
        for (int r = 0; r < 4; ++r) {
            A[r].assign(lenT, 0);
            B[r].assign(lenT, 0);
            for (int t = 0; t < lenT; ++t) {
                const int idx = (t << 2) | r;
                if (idx < na) A[r][t] = a[idx];
                if (idx < nb) B[r][t] = b[idx];
            }
        }

        std::vector<u32> r1[4], r2[4], r3[4];

        auto run_mod = [&](const NTTMod& pm, std::vector<u32> out[4]) {
            const u32 mod = pm.mod;
            for (int i = 0; i < 8; ++i) buf[i].assign((std::size_t)L, 0);

            for (int r = 0; r < 4; ++r) {
                for (int t = 0; t < lenT; ++t) {
                    buf[r][t] = (u32)(A[r][t] % (int)mod);
                    buf[4 + r][t] = (u32)(B[r][t] % (int)mod);
                }
            }

            std::vector<u32*> arrs;
            arrs.reserve(8);
            for (int i = 0; i < 8; ++i) arrs.push_back(buf[i].data());
            ntt_many(arrs, L, false, pm);

            const u32 omega = mod_pow_u32(pm.primitive_root, (mod - 1) / (u32)L, mod);
            u64 tpow = 1;

            for (int i = 0; i < L; ++i) {
                const u64 a0 = buf[0][i], a1 = buf[1][i], a2 = buf[2][i], a3 = buf[3][i];
                const u64 b0 = buf[4][i], b1 = buf[5][i], b2 = buf[6][i], b3 = buf[7][i];

                u64 c0 = a0 * b0;
                u64 c1 = a0 * b1 + a1 * b0;
                u64 c2 = a0 * b2 + a1 * b1 + a2 * b0;
                u64 c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
                u64 c4 = a1 * b3 + a2 * b2 + a3 * b1;
                u64 c5 = a2 * b3 + a3 * b2;
                u64 c6 = a3 * b3;

                c0 %= mod;
                c1 %= mod;
                c2 %= mod;
                c3 %= mod;
                c4 %= mod;
                c5 %= mod;
                c6 %= mod;

                const u64 t = tpow % mod;
                const u32 r0 = (u32)((c0 + c4 * t) % mod);
                const u32 r1v = (u32)((c1 + c5 * t) % mod);
                const u32 r2v = (u32)((c2 + c6 * t) % mod);
                const u32 r3v = (u32)c3;

                buf[0][i] = r0;
                buf[1][i] = r1v;
                buf[2][i] = r2v;
                buf[3][i] = r3v;

                tpow = (tpow * omega) % mod;
            }

            arrs.clear();
            for (int i = 0; i < 4; ++i) arrs.push_back(buf[i].data());
            ntt_many(arrs, L, true, pm);

            for (int r = 0; r < 4; ++r) {
                out[r].assign((std::size_t)lenT, 0);
                std::copy(buf[r].begin(), buf[r].begin() + lenT, out[r].begin());
            }
        };

        run_mod(P1, r1);
        run_mod(P2, r2);
        run_mod(P3, r3);

        std::vector<int> out(n_trunc, 0);
        for (int t = 0; t < lenT; ++t) {
            for (int r = 0; r < 4; ++r) {
                const int idx = (t << 2) | r;
                if (idx >= n_trunc) continue;
                out[idx] = crt3_to_mod(r1[r][t], r2[r][t], r3[r][t]);
            }
        }
        return out;
    }
};

static std::vector<u64> brute_small(int nmax) {
    std::vector<u64> dp((std::size_t)nmax + 1, 0);
    dp[0] = 1;
    for (int p = 1; p <= nmax; ++p) {
        if ((p & 3) == 2) continue;
        for (int s = nmax; s >= p; --s) dp[s] += dp[s - p];
    }
    return dp;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int N = 10'000'000;
    if (argc >= 2) N = std::stoi(argv[1]);
    assert(N >= 1000);

    const auto dp_small = brute_small(1000);
    assert(dp_small[1] == 1);
    assert(dp_small[2] == 0);
    assert(dp_small[3] == 1);
    assert(dp_small[6] == 1);
    assert(dp_small[10] == 3);
    assert(dp_small[100] == 37076);
    assert(dp_small[1000] == 3699177285485660336ULL);

    std::vector<int> invInt((std::size_t)N + 1, 0);
    invInt[1] = 1;
    for (int i = 2; i <= N; ++i) invInt[i] = MOD - (int)(1LL * (MOD / i) * invInt[MOD % i] % MOD);

    std::vector<int> a((std::size_t)N + 1, 0);
    for (int d = 1; d <= N; ++d) {
        if ((d & 3) == 2) continue;
        const int val = d % MOD;
        const int step = d << 1;
        for (int k = d; k <= N; k += step) {
            int x = a[k] + val;
            if (x >= MOD) x -= MOD;
            a[k] = x;
        }
        for (int k = step; k <= N; k += step) {
            int x = a[k] - val;
            if (x < 0) x += MOD;
            a[k] = x;
        }
    }

    std::vector<int> f((std::size_t)N + 1, 0);
    for (int k = 1; k <= N; ++k) f[k] = (int)(1LL * a[k] * invInt[k] % MOD);

    RingMul rm;
    std::vector<int> g(1, 1);
    std::vector<int> inv(1, 1);
    std::vector<int> tmp;
    std::vector<int> diff;

    for (int m = 1; m < N + 1; m <<= 1) {
        const int m2 = std::min(N + 1, m << 1);

        g.resize((std::size_t)m2, 0);
        inv.resize((std::size_t)m2, 0);

        tmp = rm.mul_trunc(g, inv, m2);
        for (int i = 0; i < m2; ++i) {
            int x = tmp[i] ? (MOD - tmp[i]) : 0;
            if (i == 0) x = (x + 2) % MOD;
            tmp[i] = x;
        }
        inv = rm.mul_trunc(inv, tmp, m2); // inv = inverse of current g (old) up to m2

        tmp.assign((std::size_t)m2 - 1, 0);
        for (int i = 1; i < m2; ++i) tmp[i - 1] = (int)(1LL * i * g[i] % MOD);
        tmp = rm.mul_trunc(tmp, inv, m2 - 1); // tmp = (g'/g) up to degree m2-2

        diff.assign((std::size_t)m2, 0);
        diff[0] = 1;
        for (int i = 1; i < m2; ++i) {
            const int log_i = (int)(1LL * tmp[i - 1] * invInt[i] % MOD);
            int x = f[i] - log_i;
            if (x < 0) x += MOD;
            diff[i] = x;
        }

        diff = rm.mul_trunc(g, diff, m2); // g_new
        g.swap(diff);

        tmp = rm.mul_trunc(g, inv, m2);
        for (int i = 0; i < m2; ++i) {
            int x = tmp[i] ? (MOD - tmp[i]) : 0;
            if (i == 0) x = (x + 2) % MOD;
            tmp[i] = x;
        }
        inv = rm.mul_trunc(inv, tmp, m2);
    }

    assert((int)g.size() == N + 1);
    assert(g[1] == 1);
    assert(g[2] == 0);
    assert(g[3] == 1);
    assert(g[6] == 1);
    assert(g[10] == 3);
    assert(g[100] == 37076);
    assert(g[1000] == (int)(dp_small[1000] % MOD));

    int ans = 0;
    for (int i = 1; i <= N; ++i) {
        ans += g[i];
        if (ans >= MOD) ans -= MOD;
    }
    std::cout << ans << "\n";
    return 0;
}
