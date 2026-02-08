#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int sum = 1000;
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

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--sum=", options.sum)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.sum > 0;
}

i64 solve(const int target_sum) {
    for (int a = 1; a < target_sum; ++a) {
        for (int b = a + 1; b < target_sum; ++b) {
            const int c = target_sum - a - b;
            if (c <= b) {
                continue;
            }
            if (a * a + b * b == c * c) {
                return static_cast<i64>(a) * static_cast<i64>(b) * static_cast<i64>(c);
            }
        }
    }

    return 0LL;
}

bool run_checkpoints() {
    if (solve(12) != 60LL) {
        std::cerr << "Checkpoint failed for sum=12" << '\n';
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

    std::cout << solve(options.sum) << '\n';
    return 0;
}
