#include <cstdint>
#include <iostream>
#include <set>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    i64 limit = 1000000000LL;
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 16;
}

i64 solve(const i64 limit) {
    i64 sum = 0;

    i64 x = 2;
    i64 y = 1;

    while (true) {
        bool any_under_limit = false;

        for (const int sign : {1, -1}) {
            const i64 num = 2 * x + sign;
            if (num % 3 != 0) {
                continue;
            }
            const i64 a = num / 3;
            if (a <= 1) {
                continue;
            }
            const i64 perimeter = 3 * a + sign;
            if (perimeter <= limit) {
                sum += perimeter;
                any_under_limit = true;
            }
        }

        const i64 next_x = 2 * x + 3 * y;
        const i64 next_y = x + 2 * y;
        x = next_x;
        y = next_y;

        if (!any_under_limit && (2 * x - 1) / 3 > limit) {
            break;
        }
    }

    return sum;
}

bool run_checkpoints() {
    if (solve(1000) != 984LL) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
