#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int power = 5;
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

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--power=", options.power)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.power >= 2;
}

i64 int_pow(const int base, const int exp) {
    i64 result = 1;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

i64 solve(const int power) {
    std::array<i64, 10> digit_power{};
    for (int d = 0; d <= 9; ++d) {
        digit_power[static_cast<std::size_t>(d)] = int_pow(d, power);
    }

    const i64 max_digit_term = digit_power[9];

    int digits = 1;
    i64 pow10 = 1;
    while (pow10 <= static_cast<i64>(digits) * max_digit_term) {
        ++digits;
        pow10 *= 10;
    }

    const i64 upper = static_cast<i64>(digits - 1) * max_digit_term;

    i64 total = 0;
    for (i64 n = 2; n <= upper; ++n) {
        i64 x = n;
        i64 sum = 0;
        while (x > 0) {
            sum += digit_power[static_cast<std::size_t>(x % 10)];
            x /= 10;
        }
        if (sum == n) {
            total += n;
        }
    }

    return total;
}

bool run_checkpoints() {
    if (solve(4) != 19316LL) {
        std::cerr << "Checkpoint failed for power=4" << '\n';
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

    std::cout << solve(options.power) << '\n';
    return 0;
}
