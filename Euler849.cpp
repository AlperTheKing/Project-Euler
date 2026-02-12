#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;

static constexpr int kMod = 1'000'000'007;

static inline int idx(int v, int s, int width) {
    return v * width + s;
}

static int F_value(int n) {
    const int S = 2 * n * (n - 1);
    const int V = 4 * (n - 1);
    const int W = S + 1;
    const int SZ = (V + 1) * W;

    std::vector<int> dp(SZ, 0), ndp(SZ, 0);
    for (int v = 0; v <= V; ++v) dp[idx(v, v, W)] = 1;

    for (int i = 2; i <= n; ++i) {
        std::fill(ndp.begin(), ndp.end(), 0);
        const int low = 2 * i * (i - 1);

        for (int t = 0; t <= S; ++t) {
            int run = 0;
            for (int v = 0; v <= V; ++v) {
                run += dp[idx(v, t, W)];
                if (run >= kMod) run -= kMod;
                int s = t + v;
                if (s <= S && s >= low) ndp[idx(v, s, W)] = run;
            }
        }
        dp.swap(ndp);
    }

    int ans = 0;
    for (int v = 0; v <= V; ++v) {
        ans += dp[idx(v, S, W)];
        if (ans >= kMod) ans -= kMod;
    }
    return ans;
}

int main() {
    if (F_value(2) != 3) return 1;
    if (F_value(7) != 32923) return 1;
    std::cout << F_value(100) << '\n';
    return 0;
}
