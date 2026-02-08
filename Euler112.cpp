#include <iostream>
#include <string>

namespace {

struct Options {
    int numerator = 99;
    int denominator = 100;
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
        if (parse_int_after_prefix(arg, "--numerator=", options.numerator) ||
            parse_int_after_prefix(arg, "--denominator=", options.denominator)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.numerator >= 1 && options.denominator >= options.numerator;
}

bool is_bouncy(int n) {
    bool increasing = false;
    bool decreasing = false;

    int prev = n % 10;
    n /= 10;

    while (n > 0) {
        const int d = n % 10;
        if (d < prev) {
            increasing = true;
        } else if (d > prev) {
            decreasing = true;
        }
        if (increasing && decreasing) {
            return true;
        }
        prev = d;
        n /= 10;
    }

    return false;
}

int solve(const int numerator, const int denominator) {
    int bouncy = 0;

    for (int n = 1;; ++n) {
        if (is_bouncy(n)) {
            ++bouncy;
        }
        if (bouncy * denominator == n * numerator) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (solve(50, 100) != 538) {
        std::cerr << "Checkpoint failed for 50 percent threshold" << '\n';
        return false;
    }
    if (solve(90, 100) != 21780) {
        std::cerr << "Checkpoint failed for 90 percent threshold" << '\n';
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

    std::cout << solve(options.numerator, options.denominator) << '\n';
    return 0;
}
