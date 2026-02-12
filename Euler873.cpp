#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <vector>

using i64 = long long;
using u64 = std::uint64_t;

static constexpr u64 MOD = 1'000'000'007ULL;

static u64 mod_pow(u64 a, u64 e) {
    u64 r = 1;
    while (e > 0) {
        if (e & 1ULL) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1ULL;
    }
    return r;
}

static u64 choose_u64(int n, int k) {
    if (k < 0 || k > n) return 0;
    k = std::min(k, n - k);
    unsigned __int128 res = 1;
    for (int i = 1; i <= k; ++i) {
        res = (res * static_cast<unsigned __int128>(n - k + i)) / static_cast<unsigned __int128>(i);
    }
    return static_cast<u64>(res);
}

static bool valid_word(const std::vector<char>& w) {
    int n = static_cast<int>(w.size());
    std::vector<int> pref(n + 1, 0);
    std::vector<int> pos_a;
    std::vector<int> pos_b;
    pos_a.reserve(n);
    pos_b.reserve(n);

    for (int i = 0; i < n; ++i) {
        pref[i + 1] = pref[i] + (w[i] == 'C' ? 1 : 0);
        if (w[i] == 'A') pos_a.push_back(i);
        if (w[i] == 'B') pos_b.push_back(i);
    }

    for (int a : pos_a) {
        for (int b : pos_b) {
            int l = std::min(a, b);
            int r = std::max(a, b);
            int cs = pref[r] - pref[l + 1];
            if (cs < 2) return false;
        }
    }
    return true;
}

static u64 brute_count(int p, int q, int r) {
    std::vector<char> w;
    w.reserve(p + q + r);
    u64 cnt = 0;

    std::function<void(int, int, int)> dfs = [&](int a, int b, int c) {
        if (a == 0 && b == 0 && c == 0) {
            if (valid_word(w)) ++cnt;
            return;
        }
        if (a > 0) {
            w.push_back('A');
            dfs(a - 1, b, c);
            w.pop_back();
        }
        if (b > 0) {
            w.push_back('B');
            dfs(a, b - 1, c);
            w.pop_back();
        }
        if (c > 0) {
            w.push_back('C');
            dfs(a, b, c - 1);
            w.pop_back();
        }
    };

    dfs(p, q, r);
    return cnt;
}

static u64 words_exact(int p, int q, int r) {
    int m = p + q;
    if (p == 0 || q == 0) return choose_u64(m + r, m);

    u64 ans = 0;
    for (int t = 1; t <= m - 1 && 2 * t <= r; ++t) {
        u64 ways_ab = 0;
        if (t & 1) {
            int u = t / 2;
            ways_ab = 2ULL * choose_u64(p - 1, u) * choose_u64(q - 1, u);
        } else {
            int u = t / 2;
            ways_ab = choose_u64(p - 1, u) * choose_u64(q - 1, u - 1)
                    + choose_u64(p - 1, u - 1) * choose_u64(q - 1, u);
        }
        if (ways_ab == 0) continue;
        ans += ways_ab * choose_u64(r - 2 * t + m, m);
    }
    return ans;
}

static u64 words_mod(i64 p, i64 q, i64 r) {
    if (p > q) std::swap(p, q);
    i64 m = p + q;

    std::vector<u64> inv(static_cast<std::size_t>(m + 1), 0);
    if (m >= 1) inv[1] = 1;
    for (i64 i = 2; i <= m; ++i) {
        inv[static_cast<std::size_t>(i)] = MOD - ((MOD / static_cast<u64>(i)) * inv[static_cast<std::size_t>(MOD % static_cast<u64>(i))]) % MOD;
    }

    u64 comb0 = 1;
    for (i64 i = 1; i <= m; ++i) {
        comb0 = (comb0 * ((r + i) % MOD)) % MOD;
        comb0 = (comb0 * inv[static_cast<std::size_t>(i)]) % MOD;
    }

    if (p == 0 || q == 0) return comb0;

    i64 max_trans = (p == q) ? (2 * p - 1) : (2 * p);
    i64 tmax = std::min({m - 1, r / 2, max_trans});
    if (tmax <= 0) return 0;

    i64 umax = tmax / 2;
    std::vector<u64> cp(static_cast<std::size_t>(umax + 1), 0);
    std::vector<u64> cq(static_cast<std::size_t>(umax + 1), 0);
    cp[0] = 1;
    cq[0] = 1;

    for (i64 u = 1; u <= umax; ++u) {
        if (u <= p - 1) {
            cp[static_cast<std::size_t>(u)] = cp[static_cast<std::size_t>(u - 1)] * ((p - u) % MOD) % MOD * inv[static_cast<std::size_t>(u)] % MOD;
        }
        if (u <= q - 1) {
            cq[static_cast<std::size_t>(u)] = cq[static_cast<std::size_t>(u - 1)] * ((q - u) % MOD) % MOD * inv[static_cast<std::size_t>(u)] % MOD;
        }
    }

    i64 n0 = r + m;
    i64 n_end = n0 - 2 * tmax;
    i64 low = n_end + 1;
    i64 high = n0;
    i64 len = high - low + 1;

    std::vector<u64> pref(static_cast<std::size_t>(len + 1), 1);
    for (i64 i = 0; i < len; ++i) {
        pref[static_cast<std::size_t>(i + 1)] = pref[static_cast<std::size_t>(i)] * ((low + i) % MOD) % MOD;
    }

    u64 inv_total = mod_pow(pref[static_cast<std::size_t>(len)], MOD - 2);
    std::vector<u64> inv_range(static_cast<std::size_t>(len), 0);
    for (i64 i = len - 1; i >= 0; --i) {
        u64 x = static_cast<u64>(low + i);
        inv_range[static_cast<std::size_t>(i)] = pref[static_cast<std::size_t>(i)] * inv_total % MOD;
        inv_total = inv_total * (x % MOD) % MOD;
    }

    auto inv_large = [&](i64 x) -> u64 {
        return inv_range[static_cast<std::size_t>(x - low)];
    };

    u64 ans = 0;
    u64 comb_t = comb0;
    i64 n_cur = n0;

    for (i64 t = 1; t <= tmax; ++t) {
        i64 n_prev = n_cur;
        comb_t = comb_t * ((n_prev - m) % MOD) % MOD;
        comb_t = comb_t * ((n_prev - m - 1) % MOD) % MOD;
        comb_t = comb_t * inv_large(n_prev) % MOD;
        comb_t = comb_t * inv_large(n_prev - 1) % MOD;
        n_cur = n_prev - 2;

        u64 ways_ab = 0;
        if (t & 1) {
            i64 u = t / 2;
            ways_ab = 2ULL * cp[static_cast<std::size_t>(u)] % MOD * cq[static_cast<std::size_t>(u)] % MOD;
        } else {
            i64 u = t / 2;
            ways_ab = (cp[static_cast<std::size_t>(u)] * cq[static_cast<std::size_t>(u - 1)]
                     + cp[static_cast<std::size_t>(u - 1)] * cq[static_cast<std::size_t>(u)]) % MOD;
        }

        ans += ways_ab * comb_t % MOD;
        if (ans >= MOD) ans -= MOD;
    }

    return ans;
}

int main() {
    assert(words_exact(2, 2, 4) == 32ULL);
    assert(words_exact(4, 4, 44) == 13'908'607'644ULL);

    for (int p = 0; p <= 3; ++p) {
        for (int q = 0; q <= 3; ++q) {
            for (int r = 0; r <= 5; ++r) {
                if (p + q + r > 9) continue;
                u64 brute = brute_count(p, q, r);
                u64 exact = words_exact(p, q, r);
                assert(brute == exact);
                assert(words_mod(p, q, r) == (exact % MOD));
            }
        }
    }

    std::cout << words_mod(1'000'000, 10'000'000, 100'000'000) << '\n';
    return 0;
}
