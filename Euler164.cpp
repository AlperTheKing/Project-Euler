#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int length = 20;
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

    return options.length >= 1;
}

u64 solve(const int length) {
    if (length == 1) {
        return 9;
    }

    std::vector<std::vector<u64>> dp(10, std::vector<u64>(10, 0));
    for (int d1 = 1; d1 <= 9; ++d1) {
        for (int d2 = 0; d2 <= 9; ++d2) {
            dp[static_cast<std::size_t>(d1)][static_cast<std::size_t>(d2)] += 1ULL;
        }
    }

    if (length == 2) {
        u64 sum = 0;
        for (int a = 0; a <= 9; ++a) {
            for (int b = 0; b <= 9; ++b) {
                sum += dp[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
            }
        }
        return sum;
    }

    for (int pos = 3; pos <= length; ++pos) {
        std::vector<std::vector<u64>> next(10, std::vector<u64>(10, 0));

        for (int a = 0; a <= 9; ++a) {
            for (int b = 0; b <= 9; ++b) {
                const u64 ways = dp[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
                if (ways == 0) {
                    continue;
                }
                for (int c = 0; c <= 9; ++c) {
                    if (a + b + c > 9) {
                        continue;
                    }
                    next[static_cast<std::size_t>(b)][static_cast<std::size_t>(c)] += ways;
                }
            }
        }

        dp = std::move(next);
    }

    u64 sum = 0;
    for (int a = 0; a <= 9; ++a) {
        for (int b = 0; b <= 9; ++b) {
            sum += dp[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(3) != 165ULL) {
        std::cerr << "Checkpoint failed for length=3" << '\n';
        return false;
    }
    if (solve(4) != 990ULL) {
        std::cerr << "Checkpoint failed for length=4" << '\n';
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
