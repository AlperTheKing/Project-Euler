#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;

static constexpr int MOD = 1'000'000'007;
static constexpr int X = 256;

static inline int mod_add(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static inline int mod_sub(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

static std::vector<int> partitions_mod(int n) {
    std::vector<int> dp(n + 1, 0);
    dp[0] = 1;
    for (int s = 1; s <= n; ++s) {
        for (int sum = s; sum <= n; ++sum) dp[sum] = mod_add(dp[sum], dp[sum - s]);
    }
    return dp;
}

static int losing_partitions_xor0(int n, const std::vector<int> &sg) {
    std::vector<std::array<int, X>> dp(n + 1);
    for (auto &row : dp) row.fill(0);
    dp[0][0] = 1;
    for (int s = 1; s <= n; ++s) {
        const int g = sg[s];
        for (int sum = s; sum <= n; ++sum) {
            auto &cur = dp[sum];
            const auto &prev = dp[sum - s];
            for (int x = 0; x < X; ++x) cur[x ^ g] = mod_add(cur[x ^ g], prev[x]);
        }
    }
    return dp[n][0];
}

static std::vector<int> sg_k2(int n) {
    std::vector<int> sg(n + 1, 0);
    for (int s = 2; s <= n; ++s) sg[s] = (s % 2 == 0 ? 1 : 0);
    return sg;
}

static std::vector<int> sg_k3(int n) {
    std::vector<int> sg(n + 1, 0);
    sg[1] = 0;
    for (int s = 2; s <= n; ++s) {
        std::array<char, X> seen{};

        for (int a = 1; a <= s - 1; ++a) seen[sg[a] ^ sg[s - a]] = 1;
        for (int a = 1; a <= s - 2; ++a) {
            for (int b = 1; b <= s - a - 1; ++b) {
                const int c = s - a - b;
                seen[sg[a] ^ sg[b] ^ sg[c]] = 1;
            }
        }

        int g = 0;
        while (g < X && seen[g]) ++g;
        sg[s] = g;
    }
    return sg;
}

static std::vector<int> sg_k4(int n) {
    std::vector<int> sg(n + 1, 0);
    for (int s = 2; s <= n; ++s) sg[s] = s - 1;
    return sg;
}

static int f(int n, const std::vector<int> &sg, int partitions_n) {
    const int losing = losing_partitions_xor0(n, sg);
    return mod_sub(partitions_n, losing);
}

static int g(int n) {
    if (n < 2) return 0;

    const auto part = partitions_mod(n);
    const int Pn = part[n];

    const auto sg2 = sg_k2(n);
    const auto sg3 = sg_k3(n);
    const auto sg4 = sg_k4(n);

    const int f2 = f(n, sg2, Pn);
    int ans = f2;
    if (n >= 3) ans = mod_add(ans, f(n, sg3, Pn));
    if (n >= 4) ans = mod_add(ans, (int)((i64)(n - 3) * f(n, sg4, Pn) % MOD));
    return ans;
}

int main() {
    {
        const int n = 5;
        const int Pn = partitions_mod(n)[n];
        assert(f(n, sg_k2(n), Pn) == 3);
        assert(f(n, sg_k3(n), Pn) == 5);
    }
    assert(g(7) == 66);
    assert(g(10) == 291);

    std::cout << g(200) << "\n";
    return 0;
}

