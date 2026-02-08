#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

bool is_prime(const std::int64_t n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (std::int64_t p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

std::int64_t solve() {
    for (int n = 9; n >= 1; --n) {
        const int digit_sum = n * (n + 1) / 2;
        if (digit_sum % 3 == 0) {
            continue;
        }

        std::string digits;
        digits.reserve(static_cast<std::size_t>(n));
        for (int d = n; d >= 1; --d) {
            digits.push_back(static_cast<char>('0' + d));
        }

        do {
            const std::int64_t value = std::stoll(digits);
            if (is_prime(value)) {
                return value;
            }
        } while (std::prev_permutation(digits.begin(), digits.end()));
    }

    return 0;
}

bool run_checkpoints() {
    if (!is_prime(2143)) {
        std::cerr << "Checkpoint failed for prime 2143" << '\n';
        return false;
    }
    if (is_prime(2145)) {
        std::cerr << "Checkpoint failed for composite 2145" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
