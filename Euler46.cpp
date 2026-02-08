#include <cmath>
#include <iostream>
#include <string>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

bool is_prime(const int n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (int p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

bool can_be_written(const int n) {
    for (int s = 1; 2 * s * s < n; ++s) {
        const int candidate = n - 2 * s * s;
        if (is_prime(candidate)) {
            return true;
        }
    }
    return false;
}

int solve() {
    for (int n = 9;; n += 2) {
        if (is_prime(n)) {
            continue;
        }
        if (!can_be_written(n)) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (!can_be_written(33)) {
        std::cerr << "Checkpoint failed for 33" << '\n';
        return false;
    }
    if (!can_be_written(45)) {
        std::cerr << "Checkpoint failed for 45" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
