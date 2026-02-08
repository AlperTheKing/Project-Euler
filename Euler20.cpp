#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int n = 100;
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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.n >= 0;
}

int solve(const int n) {
    std::vector<int> digits(1, 1);  // little-endian decimal

    for (int value = 2; value <= n; ++value) {
        int carry = 0;
        for (int& d : digits) {
            const int prod = d * value + carry;
            d = prod % 10;
            carry = prod / 10;
        }
        while (carry > 0) {
            digits.push_back(carry % 10);
            carry /= 10;
        }
    }

    int sum = 0;
    for (const int d : digits) {
        sum += d;
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(10) != 27) {
        std::cerr << "Checkpoint failed for n=10" << '\n';
        return false;
    }
    if (solve(0) != 1) {
        std::cerr << "Checkpoint failed for n=0" << '\n';
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
