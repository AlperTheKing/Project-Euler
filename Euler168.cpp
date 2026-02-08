#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int max_digits = 100;
    int modulo = 100000;
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
        if (parse_int_after_prefix(arg, "--max-digits=", options.max_digits) ||
            parse_int_after_prefix(arg, "--modulo=", options.modulo)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_digits >= 2 && options.modulo >= 1;
}

u64 search_number_mod(const int digits, const int multiplier, const int last_digit, const int modulo) {
    u64 value_mod = static_cast<u64>(last_digit % modulo);
    int current = last_digit;
    int carry = 0;
    int place_mod = 10 % modulo;
    for (int pos = 1; pos < digits; ++pos) {
        const int next = multiplier * current + carry;
        carry = next / 10;
        current = next % 10;
        value_mod = (value_mod + static_cast<u64>(current) * static_cast<u64>(place_mod)) %
                    static_cast<u64>(modulo);
        place_mod = static_cast<int>((static_cast<u64>(place_mod) * 10ULL) %
                                     static_cast<u64>(modulo));
    }

    const int leading = multiplier * current + carry;
    if (current == 0 || leading != last_digit) {
        return 0;
    }
    return value_mod;
}

u64 solve(const int max_digits, const int modulo) {
    u64 sum = 0;
    for (int digits = 2; digits <= max_digits; ++digits) {
        for (int multiplier = 1; multiplier <= 9; ++multiplier) {
            for (int last_digit = 1; last_digit <= 9; ++last_digit) {
                sum += search_number_mod(digits, multiplier, last_digit, modulo);
            }
        }
    }
    return sum % static_cast<u64>(modulo);
}

bool has_property(const int n) {
    const std::string s = std::to_string(n);
    if (s.size() <= 1) {
        return false;
    }
    if (s.back() == '0') {
        return false;
    }
    const std::string rotated = s.back() + s.substr(0, s.size() - 1);
    const int rotated_value = std::stoi(rotated);
    return (rotated_value % n) == 0;
}

u64 brute_small(const int max_digits, const int modulo) {
    int upper = 1;
    for (int i = 0; i < max_digits; ++i) {
        upper *= 10;
    }
    u64 sum = 0;
    for (int n = 11; n < upper; ++n) {
        if (has_property(n)) {
            sum = (sum + static_cast<u64>(n % modulo)) % static_cast<u64>(modulo);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (search_number_mod(6, 5, 7, 1000000000) != 142857ULL) {
        std::cerr << "Checkpoint failed for 142857 sample" << '\n';
        return false;
    }
    if (solve(6, 100000) != brute_small(6, 100000)) {
        std::cerr << "Checkpoint failed for brute cross-check at max-digits=6" << '\n';
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
    std::cout << solve(options.max_digits, options.modulo) << '\n';
    return 0;
}
