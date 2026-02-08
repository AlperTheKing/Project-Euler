#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int max_digits = 40;
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
        if (parse_int_after_prefix(arg, "--max-digits=", options.max_digits)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_digits >= 1;
}

i64 solve(const int max_digits) {
    static constexpr int kMaskAll = (1 << 10) - 1;
    std::vector<std::vector<std::array<i64, 10>>> dp(
        static_cast<std::size_t>(max_digits + 1),
        std::vector<std::array<i64, 10>>(1 << 10));

    for (int first = 1; first <= 9; ++first) {
        dp[1][1 << first][static_cast<std::size_t>(first)] = 1;
    }

    for (int len = 1; len < max_digits; ++len) {
        for (int mask = 0; mask <= kMaskAll; ++mask) {
            for (int d = 0; d <= 9; ++d) {
                const i64 ways = dp[static_cast<std::size_t>(len)][static_cast<std::size_t>(mask)][
                    static_cast<std::size_t>(d)];
                if (ways == 0) {
                    continue;
                }
                if (d > 0) {
                    const int nmask = mask | (1 << (d - 1));
                    dp[static_cast<std::size_t>(len + 1)][static_cast<std::size_t>(nmask)][
                        static_cast<std::size_t>(d - 1)] += ways;
                }
                if (d < 9) {
                    const int nmask = mask | (1 << (d + 1));
                    dp[static_cast<std::size_t>(len + 1)][static_cast<std::size_t>(nmask)][
                        static_cast<std::size_t>(d + 1)] += ways;
                }
            }
        }
    }

    i64 answer = 0;
    for (int len = 1; len <= max_digits; ++len) {
        for (int d = 0; d <= 9; ++d) {
            answer += dp[static_cast<std::size_t>(len)][kMaskAll][static_cast<std::size_t>(d)];
        }
    }
    return answer;
}

i64 brute_small(const int max_digits) {
    i64 answer = 0;
    std::function<void(int, int, int)> dfs = [&](int len, int last, int mask) {
        if (len > max_digits) {
            return;
        }
        if (mask == ((1 << 10) - 1)) {
            ++answer;
        }
        if (last > 0) {
            dfs(len + 1, last - 1, mask | (1 << (last - 1)));
        }
        if (last < 9) {
            dfs(len + 1, last + 1, mask | (1 << (last + 1)));
        }
    };
    for (int first = 1; first <= 9; ++first) {
        dfs(1, first, 1 << first);
    }
    return answer;
}

bool run_checkpoints() {
    if (solve(12) != brute_small(12)) {
        std::cerr << "Checkpoint failed for max-digits=12" << '\n';
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
    std::cout << solve(options.max_digits) << '\n';
    return 0;
}
