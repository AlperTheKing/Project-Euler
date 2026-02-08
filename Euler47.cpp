#include <iostream>
#include <string>

namespace {

struct Options {
    int consecutive = 4;
    int factors = 4;
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
        if (parse_int_after_prefix(arg, "--consecutive=", options.consecutive)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--factors=", options.factors)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.consecutive >= 1 && options.factors >= 1;
}

int distinct_prime_factor_count(int n) {
    int count = 0;

    if ((n % 2) == 0) {
        ++count;
        while ((n % 2) == 0) {
            n /= 2;
        }
    }

    for (int p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            ++count;
            while ((n % p) == 0) {
                n /= p;
            }
        }
    }

    if (n > 1) {
        ++count;
    }

    return count;
}

int solve(const int consecutive, const int factors) {
    int run = 0;
    int run_start = 2;

    for (int n = 2;; ++n) {
        if (distinct_prime_factor_count(n) == factors) {
            if (run == 0) {
                run_start = n;
            }
            ++run;
            if (run == consecutive) {
                return run_start;
            }
        } else {
            run = 0;
        }
    }
}

bool run_checkpoints() {
    if (distinct_prime_factor_count(14) != 2) {
        std::cerr << "Checkpoint failed for factor count of 14" << '\n';
        return false;
    }
    if (solve(2, 2) != 14) {
        std::cerr << "Checkpoint failed for consecutive=2 factors=2" << '\n';
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

    std::cout << solve(options.consecutive, options.factors) << '\n';
    return 0;
}
