#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int black = 60;
    int white = 40;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }

    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--black=", options.black) ||
            parse_int_after_prefix(arg, "--white=", options.white)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.black >= 0 && options.white >= 0;
}

std::vector<std::pair<int, int>> build_group_types(const int black, const int white) {
    std::vector<std::pair<int, int>> types;
    types.reserve(static_cast<std::size_t>((black + 1) * (white + 1) - 1));
    for (int b = 0; b <= black; ++b) {
        for (int w = 0; w <= white; ++w) {
            if (b == 0 && w == 0) {
                continue;
            }
            types.emplace_back(b, w);
        }
    }
    return types;
}

i64 solve(const int black, const int white) {
    const std::vector<std::pair<int, int>> types = build_group_types(black, white);

    std::vector<std::vector<i64>> dp(
        static_cast<std::size_t>(black + 1),
        std::vector<i64>(static_cast<std::size_t>(white + 1), 0));
    dp[0][0] = 1;

    for (const auto& [tb, tw] : types) {
        for (int b = black; b >= 0; --b) {
            for (int w = white; w >= 0; --w) {
                const i64 base = dp[static_cast<std::size_t>(b)][static_cast<std::size_t>(w)];
                if (base == 0) {
                    continue;
                }

                int max_repeat = 0;
                if (tb == 0) {
                    max_repeat = (tw == 0) ? 0 : (white - w) / tw;
                } else if (tw == 0) {
                    max_repeat = (black - b) / tb;
                } else {
                    max_repeat = std::min((black - b) / tb, (white - w) / tw);
                }

                for (int k = 1; k <= max_repeat; ++k) {
                    dp[static_cast<std::size_t>(b + k * tb)][static_cast<std::size_t>(w + k * tw)] += base;
                }
            }
        }
    }

    return dp[static_cast<std::size_t>(black)][static_cast<std::size_t>(white)];
}

i64 brute_dfs(const int idx, const int rem_b, const int rem_w,
              const std::vector<std::pair<int, int>>& types,
              std::unordered_map<std::uint64_t, i64>& memo) {
    if (rem_b < 0 || rem_w < 0) {
        return 0;
    }
    if (idx == static_cast<int>(types.size())) {
        return (rem_b == 0 && rem_w == 0) ? 1 : 0;
    }

    const std::uint64_t key = (static_cast<std::uint64_t>(idx) << 32) |
                              (static_cast<std::uint64_t>(rem_b) << 16) |
                              static_cast<std::uint64_t>(rem_w);
    const auto it = memo.find(key);
    if (it != memo.end()) {
        return it->second;
    }

    const int tb = types[static_cast<std::size_t>(idx)].first;
    const int tw = types[static_cast<std::size_t>(idx)].second;

    int max_repeat = 0;
    if (tb == 0) {
        max_repeat = rem_w / tw;
    } else if (tw == 0) {
        max_repeat = rem_b / tb;
    } else {
        max_repeat = std::min(rem_b / tb, rem_w / tw);
    }

    i64 ways = 0;
    for (int k = 0; k <= max_repeat; ++k) {
        ways += brute_dfs(
            idx + 1,
            rem_b - k * tb,
            rem_w - k * tw,
            types,
            memo);
    }

    memo.emplace(key, ways);
    return ways;
}

i64 brute_solve(const int black, const int white) {
    const std::vector<std::pair<int, int>> types = build_group_types(black, white);
    std::unordered_map<std::uint64_t, i64> memo;
    return brute_dfs(0, black, white, types, memo);
}

bool run_checkpoints() {
    if (solve(3, 1) != 7LL) {
        std::cerr << "Checkpoint failed for 3 black and 1 white" << '\n';
        return false;
    }
    if (solve(4, 3) != brute_solve(4, 3)) {
        std::cerr << "Checkpoint failed for brute cross-check at (4,3)" << '\n';
        return false;
    }
    if (solve(5, 2) != brute_solve(5, 2)) {
        std::cerr << "Checkpoint failed for brute cross-check at (5,2)" << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::cout << solve(options.black, options.white) << '\n';
    return 0;
}
