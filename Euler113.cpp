#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u128 = unsigned __int128;

struct Options {
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
        if (parse_int_after_prefix(arg, "--digits=", options.digits)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.digits >= 1;
}

u128 binom_u128(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    k = std::min(k, n - k);

    u128 result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * static_cast<u128>(n - k + i) / static_cast<u128>(i);
    }
    return result;
}

u128 solve_u128(const int digits) {
    const u128 increasing = binom_u128(digits + 9, 9) - 1;
    const u128 decreasing = binom_u128(digits + 10, 10) - (digits + 1);
    const u128 flat_overlap = static_cast<u128>(9 * digits);
    return increasing + decreasing - flat_overlap;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool run_checkpoints() {
    if (to_string_u128(solve_u128(6)) != "12951") {
        std::cerr << "Checkpoint failed for digits=6" << '\n';
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

    std::cout << to_string_u128(solve_u128(options.digits)) << '\n';
    return 0;
}
