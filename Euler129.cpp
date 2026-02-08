#include <iostream>
#include <numeric>
#include <string>

namespace {

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
    return options.limit >= 1;
}

int A_of_n(const int n) {
    int rem = 1 % n;
    int k = 1;
    while (rem != 0) {
        rem = (rem * 10 + 1) % n;
        ++k;
    }
    return k;
}

bool exceeds_limit(const int n, const int limit) {
    int rem = 1 % n;
    for (int k = 1; k <= limit; ++k) {
        if (rem == 0) {
            return false;
        }
        rem = (rem * 10 + 1) % n;
    }
    return true;
}

int solve(const int limit) {
    for (int n = limit + 1;; ++n) {
        if (std::gcd(n, 10) != 1) {
            continue;
        }
        if (exceeds_limit(n, limit)) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (A_of_n(7) != 6) {
        std::cerr << "Checkpoint failed for A(7)" << '\n';
        return false;
    }
    if (solve(10) != 17) {
        std::cerr << "Checkpoint failed for limit=10" << '\n';
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
