#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    i64 target = 2000000;
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

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.target >= 1;
}

i64 rectangle_count(const i64 a, const i64 b) {
    return (a * (a + 1) * b * (b + 1)) / 4;
}

i64 solve(const i64 target) {
    i64 best_area = 0;
    i64 best_diff = std::numeric_limits<i64>::max();

    for (i64 a = 1; a <= 2000; ++a) {
        for (i64 b = a; b <= 2000; ++b) {
            const i64 count = rectangle_count(a, b);
            const i64 diff = std::llabs(count - target);

            if (diff < best_diff) {
                best_diff = diff;
                best_area = a * b;
            }
            if (count > target && diff > best_diff) {
                break;
            }
        }
    }

    return best_area;
}

bool run_checkpoints() {
    if (solve(18) != 6) {
        std::cerr << "Checkpoint failed for target=18" << '\n';
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

    std::cout << solve(options.target) << '\n';
    return 0;
}
