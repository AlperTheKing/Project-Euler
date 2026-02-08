#include <boost/multiprecision/cpp_int.hpp>

#include <cmath>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_int;

struct Options {
    int limit = 100;
    int digits = 100;
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
        if (parse_int_after_prefix(arg, "--digits=", options.digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 2 && options.digits >= 1;
}

int digit_sum_sqrt_digits(const int n, const int digits) {
    cpp_int p = 0;         // current root prefix
    cpp_int remainder = 0; // current remainder

    for (int i = 0; i < digits; ++i) {
        if (i == 0) {
            remainder = n;
        } else {
            remainder *= 100;
        }

        int x = 0;
        for (int d = 9; d >= 0; --d) {
            const cpp_int probe = (20 * p + d) * d;
            if (probe <= remainder) {
                x = d;
                remainder -= probe;
                break;
            }
        }

        p = p * 10 + x;
    }

    const std::string s = p.convert_to<std::string>();
    int sum = 0;
    for (char c : s) {
        sum += c - '0';
    }
    return sum;
}

int solve(const int limit, const int digits) {
    int total = 0;
    for (int n = 1; n <= limit; ++n) {
        const int r = static_cast<int>(std::sqrt(n));
        if (r * r == n) {
            continue;
        }
        total += digit_sum_sqrt_digits(n, digits);
    }
    return total;
}

bool run_checkpoints() {
    if (digit_sum_sqrt_digits(2, 100) != 475) {
        std::cerr << "Checkpoint failed for sqrt(2) first 100 digits" << '\n';
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

    std::cout << solve(options.limit, options.digits) << '\n';
    return 0;
}
