#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i128 = __int128_t;

struct Options {
    int turns = 15;
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
        if (parse_int_after_prefix(arg, "--turns=", options.turns)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.turns >= 1;
}

i128 prize_fund(const int turns) {
    std::vector<i128> dp(static_cast<std::size_t>(turns + 1), 0);
    dp[0] = 1;

    for (int t = 1; t <= turns; ++t) {
        std::vector<i128> next(static_cast<std::size_t>(turns + 1), 0);
        for (int blue = 0; blue <= t; ++blue) {
            if (blue <= t - 1) {
                next[static_cast<std::size_t>(blue)] += dp[static_cast<std::size_t>(blue)] * t;
            }
            if (blue > 0) {
                next[static_cast<std::size_t>(blue)] += dp[static_cast<std::size_t>(blue - 1)];
            }
        }
        dp = std::move(next);
    }

    i128 denominator = 1;
    for (int t = 1; t <= turns; ++t) {
        denominator *= (t + 1);
    }

    i128 win_numerator = 0;
    for (int blue = turns / 2 + 1; blue <= turns; ++blue) {
        win_numerator += dp[static_cast<std::size_t>(blue)];
    }

    return denominator / win_numerator;
}

bool run_checkpoints() {
    if (prize_fund(4) != 10) {
        std::cerr << "Checkpoint failed for turns=4" << '\n';
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

    std::cout << static_cast<long long>(prize_fund(options.turns)) << '\n';
    return 0;
}
