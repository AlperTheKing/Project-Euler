#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

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

bool has_substring_divisibility(const std::string& s) {
    static const int divisors[7] = {2, 3, 5, 7, 11, 13, 17};

    if (s.size() != 10U) {
        return false;
    }
    if (s[0] == '0') {
        return false;
    }

    for (int i = 0; i < 7; ++i) {
        const int value = (s[static_cast<std::size_t>(i + 1)] - '0') * 100 +
                          (s[static_cast<std::size_t>(i + 2)] - '0') * 10 +
                          (s[static_cast<std::size_t>(i + 3)] - '0');
        if (value % divisors[i] != 0) {
            return false;
        }
    }

    return true;
}

i64 solve() {
    std::string digits = "0123456789";
    i64 total = 0;

    do {
        if (has_substring_divisibility(digits)) {
            total += std::stoll(digits);
        }
    } while (std::next_permutation(digits.begin(), digits.end()));

    return total;
}

bool run_checkpoints() {
    if (!has_substring_divisibility("1406357289")) {
        std::cerr << "Checkpoint failed for 1406357289" << '\n';
        return false;
    }
    if (has_substring_divisibility("1234567890")) {
        std::cerr << "Checkpoint failed for invalid example" << '\n';
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
