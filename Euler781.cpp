#include <cstdint>
#include <iostream>
#include <vector>

namespace {

constexpr int MOD = 1000000007;

int solve_dp(int n) {
    if (n & 1) return 0;
    std::vector<int> dp{0, 0, 1};
    dp.reserve(static_cast<std::size_t>(n + 2));

    for (int i = 4; i <= n; i += 2) {
        std::vector<int> next = dp;
        for (int j = static_cast<int>(next.size()) - 2; j >= 1; --j) {
            int v = next[j] + next[j + 1];
            if (v >= MOD) v -= MOD;
            next[j] = v;
        }
        next.resize(dp.size() + 2, 0);
        for (int j = 1; j < static_cast<int>(dp.size()); ++j) {
            const std::int64_t add = static_cast<std::int64_t>(j + 1) * dp[j] % MOD;
            int idx = j + 2;
            int v = next[idx] + static_cast<int>(add);
            if (v >= MOD) v -= MOD;
            next[idx] = v;
        }
        dp.swap(next);
    }

    std::int64_t sum = 0;
    for (int v : dp) {
        sum += v;
        if (sum >= MOD) sum -= MOD;
    }
    return static_cast<int>(sum % MOD);
}

}  // namespace

int main(int argc, char** argv) {
    int n = 50000;
    if (argc >= 2) n = std::stoi(argv[1]);

    if (n >= 4 && solve_dp(4) != 5) {
        std::cerr << "Validation failed for F(4)." << '\n';
        return 1;
    }
    if (n >= 8 && solve_dp(8) != 319) {
        std::cerr << "Validation failed for F(8)." << '\n';
        return 1;
    }

    std::cout << solve_dp(n) << '\n';
    return 0;
}
