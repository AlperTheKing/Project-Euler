#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int exponent = 1000;
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
        if (parse_int_after_prefix(arg, "--exponent=", options.exponent)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.exponent >= 0;
}

int solve(const int exponent) {
    std::vector<int> digits(1, 1);  // little-endian decimal

    for (int e = 0; e < exponent; ++e) {
        int carry = 0;
        for (int& d : digits) {
            const int value = 2 * d + carry;
            d = value % 10;
            carry = value / 10;
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
    if (solve(15) != 26) {
        std::cerr << "Checkpoint failed for exponent=15" << '\n';
        return false;
    }
    if (solve(0) != 1) {
        std::cerr << "Checkpoint failed for exponent=0" << '\n';
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

    std::cout << solve(options.exponent) << '\n';
    return 0;
}
