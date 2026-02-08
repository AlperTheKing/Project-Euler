#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int bits = 32;
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
        if (parse_int_after_prefix(arg, "--bits=", options.bits)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.bits >= 1 && options.bits <= 128;
}

long double solve(const int bits) {
    long double expected = 0.0L;
    for (int n = 0; n < 10000; ++n) {
        const long double bit_zero_prob = std::ldexp(1.0L, -n);  // 2^-n
        const long double all_ones_prob = std::powl(1.0L - bit_zero_prob, static_cast<long double>(bits));
        const long double term = 1.0L - all_ones_prob;
        expected += term;
        if (term < 1e-20L) {
            break;
        }
    }
    return expected;
}

bool run_checkpoints() {
    if (std::fabsl(solve(1) - 2.0L) > 1e-15L) {
        std::cerr << "Checkpoint failed for 1-bit case" << '\n';
        return false;
    }
    if (std::fabsl(solve(2) - (8.0L / 3.0L)) > 1e-15L) {
        std::cerr << "Checkpoint failed for 2-bit case" << '\n';
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
    const long double answer = solve(options.bits);
    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(answer) << '\n';
    return 0;
}
