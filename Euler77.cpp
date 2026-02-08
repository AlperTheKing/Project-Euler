#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int target = 5000;
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

int solve(const int target) {
    std::vector<int> primes;

    for (int n = 2;; ++n) {
        if (is_prime(n)) {
            primes.push_back(n);
        }

        std::vector<i64> ways(static_cast<std::size_t>(n + 1), 0LL);
        ways[0] = 1;

        for (const int p : primes) {
            if (p > n) {
                break;
            }
            for (int sum = p; sum <= n; ++sum) {
                ways[static_cast<std::size_t>(sum)] += ways[static_cast<std::size_t>(sum - p)];
            }
        }

        if (ways[static_cast<std::size_t>(n)] > target) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (solve(5) != 11) {
        std::cerr << "Checkpoint failed for target=5" << '\n';
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
