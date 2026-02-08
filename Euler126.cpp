#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int target = 1000;
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

int layer_size(const int a, const int b, const int c, const int n) {
    return 2 * (a * b + a * c + b * c) +
           4 * (n - 1) * (a + b + c) +
           4 * (n - 1) * (n - 2);
}

int solve(const int target) {
    int limit = 50000;

    for (;;) {
        std::vector<int> count(static_cast<std::size_t>(limit + 1), 0);

        for (int a = 1; layer_size(a, a, a, 1) <= limit; ++a) {
            for (int b = a; layer_size(a, b, b, 1) <= limit; ++b) {
                for (int c = b; layer_size(a, b, c, 1) <= limit; ++c) {
                    for (int n = 1;; ++n) {
                        const int cubes = layer_size(a, b, c, n);
                        if (cubes > limit) {
                            break;
                        }
                        ++count[static_cast<std::size_t>(cubes)];
                    }
                }
            }
        }

        for (int i = 1; i <= limit; ++i) {
            if (count[static_cast<std::size_t>(i)] == target) {
                return i;
            }
        }

        limit *= 2;
    }
}

bool run_checkpoints() {
    if (layer_size(1, 1, 1, 1) != 6) {
        std::cerr << "Checkpoint failed for C(1,1,1,1)" << '\n';
        return false;
    }
    if (solve(10) != 154) {
        std::cerr << "Checkpoint failed for target=10" << '\n';
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
