#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

static constexpr u32 kMod = 999'676'999U;

static std::vector<int> build_spf(int n) {
    std::vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) spf[i] = i;
    for (int i = 2; static_cast<u64>(i) * i <= static_cast<u64>(n); ++i) {
        if (spf[i] != i) continue;
        for (int j = i * i; j <= n; j += i) {
            if (spf[j] == j) spf[j] = i;
        }
    }
    return spf;
}

static std::vector<u32> build_D_mod(int n, const std::vector<int>& spf) {
    std::vector<u32> d(n + 1, 0);
    d[1] = 1;
    for (int x = 2; x <= n; ++x) {
        if (spf[x] == x) {
            d[x] = 1;
        } else {
            const int p = spf[x];
            const int q = x / p;
            d[x] = static_cast<u32>((static_cast<u64>(p) * d[q] + q) % kMod);
        }
    }
    return d;
}

static std::vector<u32> build_G_mod(int n, const std::vector<u32>& d) {
    std::vector<u32> g(n + 1, 0);
    g[0] = 1;
    for (int part = 1; part <= n; ++part) {
        const u32 w = d[part];
        if (w == 1U) {
            for (int s = part; s <= n; ++s) {
                u32 v = g[s] + g[s - part];
                if (v >= kMod) v -= kMod;
                g[s] = v;
            }
        } else {
            for (int s = part; s <= n; ++s) {
                g[s] = static_cast<u32>((g[s] + static_cast<u64>(w) * g[s - part]) % kMod);
            }
        }
    }
    return g;
}

int main() {
    const int n = 50'000;
    const std::vector<int> spf = build_spf(n);
    const std::vector<u32> d = build_D_mod(n, spf);
    const std::vector<u32> g = build_G_mod(n, d);

    assert(g[10] == 164U);
    u32 s10 = 0;
    for (int i = 1; i <= 10; ++i) {
        s10 += g[i];
        if (s10 >= kMod) s10 -= kMod;
    }
    assert(s10 == 396U);

    u32 ans = 0;
    for (int i = 1; i <= n; ++i) {
        ans += g[i];
        if (ans >= kMod) ans -= kMod;
    }
    std::cout << ans << '\n';
    return 0;
}
