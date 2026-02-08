#include <cmath>
#include <iostream>
#include <string>

namespace {

struct Options {
    int target = 1000000;
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
        if (parse_int_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 1;
}

int solve(const int target) {
    int total = 0;

    for (int m = 1;; ++m) {
        for (int s = 2; s <= 2 * m; ++s) {
            const int value = m * m + s * s;
            const int r = static_cast<int>(std::sqrt(static_cast<double>(value)));
            if (r * r != value) {
                continue;
            }

            const int lo = std::max(1, s - m);
            const int hi = std::min(m, s / 2);
            if (hi >= lo) {
                total += hi - lo + 1;
            }
        }

        if (total > target) {
            return m;
        }
    }
}

bool run_checkpoints() {
    if (solve(2000) != 100) {
        std::cerr << "Checkpoint failed for target=2000" << '\n';
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
