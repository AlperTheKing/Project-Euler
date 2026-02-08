#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int count = 30;
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
        if (parse_int_after_prefix(arg, "--count=", options.count)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.count >= 1;
}

i64 solve(const int count) {
    std::vector<std::pair<i64, i64>> seeds = {
        {7, 1}, {8, 2}, {13, 5}, {17, 7}, {32, 14}, {43, 19},
    };

    std::set<i64> nuggets;

    for (auto [x0, y0] : seeds) {
        i64 x = x0;
        i64 y = y0;

        for (int step = 0; step < 50; ++step) {
            if (x > 7 && (x - 7) % 5 == 0) {
                nuggets.insert((x - 7) / 5);
            }

            const i64 next_x = 9 * x + 20 * y;
            const i64 next_y = 4 * x + 9 * y;
            x = next_x;
            y = next_y;
        }
    }

    std::vector<i64> values(nuggets.begin(), nuggets.end());
    std::sort(values.begin(), values.end());

    i64 sum = 0;
    for (int i = 0; i < count; ++i) {
        sum += values[static_cast<std::size_t>(i)];
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(5) != 222) {
        std::cerr << "Checkpoint failed for first five nuggets" << '\n';
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

    std::cout << solve(options.count) << '\n';
    return 0;
}
