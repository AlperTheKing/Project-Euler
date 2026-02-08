#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    i64 limit = 10000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 1;
}

int solve(const i64 limit) {
    std::vector<int> primes;
    primes.reserve(300000);

    for (int candidate = 2;; ++candidate) {
        bool prime = true;
        for (int p : primes) {
            if (static_cast<i64>(p) * p > candidate) {
                break;
            }
            if (candidate % p == 0) {
                prime = false;
                break;
            }
        }
        if (!prime) {
            continue;
        }

        primes.push_back(candidate);
        const int n = static_cast<int>(primes.size());
        if ((n & 1) == 0) {
            continue;
        }

        const i64 p = candidate;
        const i64 remainder = 2LL * static_cast<i64>(n) * p;
        if (remainder > limit) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (solve(1000000000LL) != 7037) {
        std::cerr << "Checkpoint failed for limit=1e9" << '\n';
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
