#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

i64 G(const int n) {
    const int mod = n + 1;
    std::vector<int> pos(static_cast<std::size_t>(n + 1), 0);

    int value = 1;
    for (int i = 1; i <= n; ++i) {
        value = static_cast<int>((3LL * value) % mod);
        if (value <= 0 || value > n || pos[static_cast<std::size_t>(value)] != 0) {
            return -1;
        }
        pos[static_cast<std::size_t>(value)] = i;
    }

    for (int card = 1; card <= n; ++card) {
        if (pos[static_cast<std::size_t>(card)] == 0) {
            return -1;
        }
    }

    const i64 inf = (1LL << 60);
    std::vector<std::vector<i64>> dp(static_cast<std::size_t>(n + 1),
                                     std::vector<i64>(static_cast<std::size_t>(n + 1), 0));

    for (int len = 2; len <= n; ++len) {
        for (int l = 1; l + len - 1 <= n; ++l) {
            const int r = l + len - 1;
            i64 best = inf;
            const int pos_r = pos[static_cast<std::size_t>(r)];
            for (int m = l; m < r; ++m) {
                const i64 candidate =
                    dp[static_cast<std::size_t>(l)][static_cast<std::size_t>(m)] +
                    dp[static_cast<std::size_t>(m + 1)][static_cast<std::size_t>(r)] +
                    std::llabs(static_cast<i64>(pos[static_cast<std::size_t>(m)]) -
                               static_cast<i64>(pos_r));
                if (candidate < best) {
                    best = candidate;
                }
            }
            dp[static_cast<std::size_t>(l)][static_cast<std::size_t>(r)] = best;
        }
    }

    return dp[1][static_cast<std::size_t>(n)];
}

}  // namespace

int main() {
    assert(G(6) == 8);
    assert(G(16) == 47);

    std::cout << G(976) << '\n';
    return 0;
}
