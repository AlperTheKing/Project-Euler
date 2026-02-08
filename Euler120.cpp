#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int max_a = 1000;
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
        if (parse_int_after_prefix(arg, "--max-a=", options.max_a)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.max_a >= 3;
}

i64 r_max(const int a) {
    if ((a & 1) == 0) {
        return static_cast<i64>(a) * static_cast<i64>(a - 2);
    }
    return static_cast<i64>(a) * static_cast<i64>(a - 1);
}

i64 solve(const int max_a) {
    i64 total = 0;
    for (int a = 3; a <= max_a; ++a) {
        total += r_max(a);
    }
    return total;
}

bool run_checkpoints() {
    if (r_max(7) != 42) {
        std::cerr << "Checkpoint failed for a=7" << '\n';
        return false;
    }
    if (solve(7) != 100) {
        std::cerr << "Checkpoint failed for range 3..7" << '\n';
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

    std::cout << solve(options.max_a) << '\n';
    return 0;
}
