#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 12000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2;
}

i64 count_between(const int a, const int b, const int c, const int d, const int limit) {
    const int num = a + c;
    const int den = b + d;
    if (den > limit) {
        return 0;
    }

    return 1 + count_between(a, b, num, den, limit) + count_between(num, den, c, d, limit);
}

i64 solve(const int limit) {
    return count_between(1, 3, 1, 2, limit);
}

bool run_checkpoints() {
    if (solve(8) != 3LL) {
        std::cerr << "Checkpoint failed for limit=8" << '\n';
        return false;
    }
    if (solve(12) != 7LL) {
        std::cerr << "Checkpoint failed for limit=12" << '\n';
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
