#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 MOD = 987'898'789U;
constexpr int TARGET_N = 500;
constexpr std::array<int, 5> K_VALUES = {1, 10, 100, 1000, 10000};

u32 mod_add(const u32 a, const u32 b) {
    const u32 sum = a + b;
    return sum >= MOD ? sum - MOD : sum;
}

u32 mod_mul(const u32 a, const u32 b) {
    return static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(b)) % static_cast<u64>(MOD));
}

std::vector<std::vector<u32>> build_needed_binomials(const int n, const std::array<int, 5>& k_values) {
    int max_row = 0;
    std::vector<bool> needed(1, false);
    for (const int k : k_values) {
        max_row = std::max(max_row, k + n - 2);
    }
    needed.assign(static_cast<std::size_t>(max_row + 1), false);
    for (const int k : k_values) {
        for (int row = k; row <= k + n - 2; ++row) {
            needed[static_cast<std::size_t>(row)] = true;
        }
    }

    std::vector<std::vector<u32>> rows(static_cast<std::size_t>(max_row + 1));
    std::vector<u32> current(1, 1U);
    for (int row = 1; row <= max_row; ++row) {
        current.push_back(1U);
        for (int col = row - 1; col >= 1; --col) {
            current[static_cast<std::size_t>(col)] =
                mod_add(current[static_cast<std::size_t>(col)], current[static_cast<std::size_t>(col - 1)]);
        }
        if (needed[static_cast<std::size_t>(row)]) {
            rows[static_cast<std::size_t>(row)] = current;
        }
    }
    return rows;
}

u32 solve_endpoint(const int n, const int k, const int endpoint, const std::vector<std::vector<u32>>& binomials) {
    int right_crossings = k - (endpoint == 0 ? 1 : 0);
    u32 ways = 1U;

    for (int stone = 1; stone < n; ++stone) {
        const int row = k + stone - 1;
        if (right_crossings <= 0 || right_crossings > row + 1) {
            return 0U;
        }
        ways = mod_mul(ways, binomials[static_cast<std::size_t>(row)][static_cast<std::size_t>(right_crossings - 1)]);
        right_crossings = k + stone - right_crossings + (endpoint > stone ? 1 : 0);
    }

    return ways;
}

u32 solve(const int n, const int k, const std::vector<std::vector<u32>>& binomials) {
    u32 total = 0U;
    for (int endpoint = 0; endpoint <= n; ++endpoint) {
        total = mod_add(total, solve_endpoint(n, k, endpoint, binomials));
    }
    return total;
}

u64 brute_force(const int n, const int k) {
    const std::vector<int> target = [&]() {
        std::vector<int> counts(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            counts[static_cast<std::size_t>(i)] = k + i;
        }
        return counts;
    }();

    std::map<std::vector<int>, u64> memo;

    std::function<u64(int, std::vector<int>&)> dfs = [&](const int position, std::vector<int>& counts) -> u64 {
        std::vector<int> state;
        state.reserve(static_cast<std::size_t>(n + 1));
        state.push_back(position);
        state.insert(state.end(), counts.begin(), counts.end());

        const auto it = memo.find(state);
        if (it != memo.end()) {
            return it->second;
        }

        bool complete = true;
        for (int i = 0; i < n; ++i) {
            if (counts[static_cast<std::size_t>(i)] > target[static_cast<std::size_t>(i)]) {
                return 0ULL;
            }
            if (counts[static_cast<std::size_t>(i)] != target[static_cast<std::size_t>(i)]) {
                complete = false;
            }
        }

        if (complete) {
            const u64 total = 1ULL + (position == n - 1 ? 1ULL : 0ULL);
            memo.emplace(std::move(state), total);
            return total;
        }

        u64 total = 0ULL;
        for (const int next : {position - 1, position + 1}) {
            if (next < 0 || next > n) {
                continue;
            }
            if (next < n) {
                ++counts[static_cast<std::size_t>(next)];
                if (counts[static_cast<std::size_t>(next)] <= target[static_cast<std::size_t>(next)]) {
                    total += dfs(next, counts);
                }
                --counts[static_cast<std::size_t>(next)];
            } else {
                total += dfs(next, counts);
            }
        }

        memo.emplace(std::move(state), total);
        return total;
    };

    std::vector<int> counts(static_cast<std::size_t>(n), 0);
    counts[0] = 1;
    return dfs(0, counts);
}

void run_checkpoints(const std::vector<std::vector<u32>>& binomials) {
    for (int n = 1; n <= 6; ++n) {
        for (int k = 1; k <= 3; ++k) {
            assert(static_cast<u64>(solve(n, k, binomials)) == brute_force(n, k));
        }
    }

    assert(solve(3, 2, binomials) == 17U);
    assert(solve(6, 1, binomials) == 1320U);
    assert(solve(6, 5, binomials) == 16'793'280U);
}

u32 solve_target(const std::vector<std::vector<u32>>& binomials) {
    u32 total = 0U;
    for (const int k : K_VALUES) {
        total = mod_add(total, solve(TARGET_N, k, binomials));
    }
    return total;
}

}  // namespace

int main(int argc, char** argv) {
    bool should_run_checkpoints = true;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--skip-checkpoints") {
            should_run_checkpoints = false;
        }
    }

    const std::vector<std::vector<u32>> binomials = build_needed_binomials(TARGET_N, K_VALUES);

    if (should_run_checkpoints) {
        run_checkpoints(binomials);
    }

    std::cout << solve_target(binomials) << '\n';
    return 0;
}
