#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int days = 30;
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
        if (parse_int_after_prefix(arg, "--days=", options.days)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.days >= 1;
}

i64 solve(const int days) {
    // dp[late_used][consecutive_absent]
    i64 dp[2][3] = {};
    dp[0][0] = 1;

    for (int day = 0; day < days; ++day) {
        i64 next[2][3] = {};
        for (int late = 0; late <= 1; ++late) {
            for (int a = 0; a <= 2; ++a) {
                const i64 ways = dp[late][a];
                if (ways == 0) {
                    continue;
                }

                // On time
                next[late][0] += ways;
                // Absent
                if (a < 2) {
                    next[late][a + 1] += ways;
                }
                // Late
                if (late == 0) {
                    next[1][0] += ways;
                }
            }
        }

        for (int late = 0; late <= 1; ++late) {
            for (int a = 0; a <= 2; ++a) {
                dp[late][a] = next[late][a];
            }
        }
    }

    i64 total = 0;
    for (int late = 0; late <= 1; ++late) {
        for (int a = 0; a <= 2; ++a) {
            total += dp[late][a];
        }
    }
    return total;
}

bool run_checkpoints() {
    if (solve(4) != 43LL) {
        std::cerr << "Checkpoint failed for days=4" << '\n';
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

    std::cout << solve(options.days) << '\n';
    return 0;
}
