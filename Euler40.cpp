#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int max_power = 6;
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
    for (const char c : tail) {
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
        if (parse_int_after_prefix(arg, "--max-power=", options.max_power)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.max_power >= 0;
}

i64 solve(const int max_power) {
    std::vector<int> targets;
    targets.reserve(static_cast<std::size_t>(max_power + 1));

    int t = 1;
    for (int p = 0; p <= max_power; ++p) {
        targets.push_back(t);
        if (t > std::numeric_limits<int>::max() / 10) {
            throw std::overflow_error("target overflow");
        }
        t *= 10;
    }

    const int max_target = targets.back();

    i64 product = 1;
    int next_target_index = 0;
    int global_pos = 0;

    for (int n = 1; global_pos < max_target; ++n) {
        const std::string s = std::to_string(n);
        for (const char c : s) {
            ++global_pos;
            if (global_pos == targets[static_cast<std::size_t>(next_target_index)]) {
                product *= static_cast<i64>(c - '0');
                ++next_target_index;
                if (next_target_index >= static_cast<int>(targets.size())) {
                    return product;
                }
            }
        }
    }

    return product;
}

bool run_checkpoints() {
    if (solve(2) != 5LL) {
        std::cerr << "Checkpoint failed for max_power=2" << '\n';
        return false;
    }
    if (solve(0) != 1LL) {
        std::cerr << "Checkpoint failed for max_power=0" << '\n';
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

    try {
        std::cout << solve(options.max_power) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
