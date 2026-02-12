#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;

static constexpr i64 kMod = 998244353;

enum Dir : int { L = 0, R = 1, V = 2 };

static int swap_low2(int mask) {
    return ((mask & 1) << 1) | ((mask >> 1) & 1);
}

static std::array<std::vector<std::array<int, 3>>, 4> build_transitions(bool last_column) {
    std::array<std::vector<std::array<int, 3>>, 4> grouped;

    for (int in = 0; in < 4; ++in) {
        std::array<std::array<int, 3>, 4> cnt{};
        const int in_top = in & 1;
        const int in_bottom = (in >> 1) & 1;

        for (int top = 0; top < 3; ++top) {
            for (int bot = 0; bot < 3; ++bot) {
                int delta = 0;
                if (top == L && in_top) {
                    ++delta;
                }
                if (bot == L && in_bottom) {
                    ++delta;
                }
                if (top == V && bot == V) {
                    ++delta;
                }

                int out = 0;
                if (top == R) {
                    out |= 1;
                }
                if (bot == R) {
                    out |= 2;
                }
                if (last_column) {
                    out = swap_low2(out);
                }
                ++cnt[out][delta];
            }
        }

        for (int out = 0; out < 4; ++out) {
            for (int d = 0; d <= 2; ++d) {
                if (cnt[out][d] > 0) {
                    grouped[in].push_back({out, d, cnt[out][d]});
                }
            }
        }
    }

    return grouped;
}

static i64 S_mod(int n) {
    const int cols = 2 * n;
    const int target = n;

    const auto trans_mid = build_transitions(false);
    const auto trans_last = build_transitions(true);

    i64 total = 0;

    for (int init = 0; init < 4; ++init) {
        std::array<std::vector<i64>, 4> dp;
        std::array<std::vector<i64>, 4> next;
        for (int s = 0; s < 4; ++s) {
            dp[s].assign(target + 1, 0);
            next[s].assign(target + 1, 0);
        }
        dp[init][0] = 1;

        for (int col = 0; col < cols - 1; ++col) {
            for (int s = 0; s < 4; ++s) {
                std::fill(next[s].begin(), next[s].end(), 0);
            }

            for (int in = 0; in < 4; ++in) {
                const auto& cur = dp[in];
                for (const auto& tr : trans_mid[in]) {
                    const int out = tr[0];
                    const int delta = tr[1];
                    const int ways = tr[2];

                    for (int k = 0; k + delta <= target; ++k) {
                        if (cur[k] == 0) {
                            continue;
                        }
                        i64 add = (cur[k] * ways) % kMod;
                        i64& dst = next[out][k + delta];
                        dst += add;
                        if (dst >= kMod) {
                            dst -= kMod;
                        }
                    }
                }
            }

            dp.swap(next);
        }

        for (int s = 0; s < 4; ++s) {
            std::fill(next[s].begin(), next[s].end(), 0);
        }

        for (int in = 0; in < 4; ++in) {
            const auto& cur = dp[in];
            for (const auto& tr : trans_last[in]) {
                const int out = tr[0];
                const int delta = tr[1];
                const int ways = tr[2];

                for (int k = 0; k + delta <= target; ++k) {
                    if (cur[k] == 0) {
                        continue;
                    }
                    i64 add = (cur[k] * ways) % kMod;
                    i64& dst = next[out][k + delta];
                    dst += add;
                    if (dst >= kMod) {
                        dst -= kMod;
                    }
                }
            }
        }

        total += next[init][target];
        total %= kMod;
    }

    return total;
}

int main() {
    assert(S_mod(1) == 48);
    assert(S_mod(10) == 420121075);
    std::cout << S_mod(1000) << '\n';
    return 0;
}
