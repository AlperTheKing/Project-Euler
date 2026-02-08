#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int size = 1001;
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
        if (parse_int_after_prefix(arg, "--size=", options.size)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.size > 0 && (options.size % 2 == 1);
}

i64 solve(const int size) {
    i64 total = 1;
    for (int side = 3; side <= size; side += 2) {
        total += 4LL * side * side - 6LL * (side - 1);
    }
    return total;
}

bool run_checkpoints() {
    if (solve(5) != 101LL) {
        std::cerr << "Checkpoint failed for size=5" << '\n';
        return false;
    }
    if (solve(7) != 261LL) {
        std::cerr << "Checkpoint failed for size=7" << '\n';
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

    std::cout << solve(options.size) << '\n';
    return 0;
}
