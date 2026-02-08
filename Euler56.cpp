#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int a_max = 99;
    int b_max = 99;
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
        if (parse_int_after_prefix(arg, "--a-max=", options.a_max)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--b-max=", options.b_max)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.a_max >= 1 && options.b_max >= 1;
}

void multiply_small(std::vector<int>& digits, const int mul) {
    int carry = 0;
    for (int& d : digits) {
        const int value = d * mul + carry;
        d = value % 10;
        carry = value / 10;
    }
    while (carry > 0) {
        digits.push_back(carry % 10);
        carry /= 10;
    }
}

int digit_sum(const std::vector<int>& digits) {
    int sum = 0;
    for (const int d : digits) {
        sum += d;
    }
    return sum;
}

int solve(const int a_max, const int b_max) {
    int best = 0;

    for (int a = 1; a <= a_max; ++a) {
        std::vector<int> value(1, 1);
        for (int b = 1; b <= b_max; ++b) {
            multiply_small(value, a);
            best = std::max(best, digit_sum(value));
        }
    }

    return best;
}

bool run_checkpoints() {
    if (solve(10, 10) != 45) {
        std::cerr << "Checkpoint failed for 10x10 range" << '\n';
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

    std::cout << solve(options.a_max, options.b_max) << '\n';
    return 0;
}
