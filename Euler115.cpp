#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int min_block = 50;
    i64 target = 1000000;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
    }
    value = parsed;
    return true;
}

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
        if (parse_i64_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        if (parse_int_after_prefix(arg, "--min-block=", options.min_block)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.min_block >= 2 && options.target >= 1;
}

i64 count_ways(const int length, const int min_block) {
    std::vector<i64> dp(static_cast<std::size_t>(length + 1), 0);
    dp[0] = 1;

    for (int n = 1; n <= length; ++n) {
        i64 ways = dp[static_cast<std::size_t>(n - 1)];

        for (int block = min_block; block <= n; ++block) {
            if (block == n) {
                ways += 1;
            } else {
                ways += dp[static_cast<std::size_t>(n - block - 1)];
            }
        }

        dp[static_cast<std::size_t>(n)] = ways;
    }

    return dp[static_cast<std::size_t>(length)];
}

int solve(const int min_block, const i64 target) {
    for (int n = 0;; ++n) {
        if (count_ways(n, min_block) > target) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (solve(3, 1000000) != 30) {
        std::cerr << "Checkpoint failed for min_block=3" << '\n';
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

    std::cout << solve(options.min_block, options.target) << '\n';
    return 0;
}
