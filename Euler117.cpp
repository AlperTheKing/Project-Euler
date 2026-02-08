#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int length = 50;
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
        if (parse_int_after_prefix(arg, "--length=", options.length)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.length >= 0;
}

i64 solve(const int length) {
    std::vector<i64> dp(static_cast<std::size_t>(length + 1), 0);
    dp[0] = 1;

    for (int n = 1; n <= length; ++n) {
        i64 ways = dp[static_cast<std::size_t>(n - 1)];
        if (n >= 2) {
            ways += dp[static_cast<std::size_t>(n - 2)];
        }
        if (n >= 3) {
            ways += dp[static_cast<std::size_t>(n - 3)];
        }
        if (n >= 4) {
            ways += dp[static_cast<std::size_t>(n - 4)];
        }
        dp[static_cast<std::size_t>(n)] = ways;
    }

    return dp[static_cast<std::size_t>(length)];
}

bool run_checkpoints() {
    if (solve(5) != 15) {
        std::cerr << "Checkpoint failed for length=5" << '\n';
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

    std::cout << solve(options.length) << '\n';
    return 0;
}
