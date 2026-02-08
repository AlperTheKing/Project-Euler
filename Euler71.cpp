#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 1000000;
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

int solve(const int limit) {
    i64 best_n = 0;
    i64 best_d = 1;

    for (int d = 2; d <= limit; ++d) {
        const i64 n = (3LL * d - 1) / 7LL;
        if (std::gcd(static_cast<i64>(d), n) != 1LL) {
            continue;
        }
        if (n * best_d > best_n * d) {
            best_n = n;
            best_d = d;
        }
    }

    return static_cast<int>(best_n);
}

bool run_checkpoints() {
    if (solve(8) != 2) {
        std::cerr << "Checkpoint failed for limit=8" << '\n';
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
