#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int exponent = 30;  // n <= 2^exponent
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
        if (parse_int_after_prefix(arg, "--exponent=", options.exponent)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.exponent >= 1 && options.exponent <= 60;
}

u64 solve(const int exponent) {
    const int bits = exponent + 1;
    std::vector<u64> fib(static_cast<std::size_t>(bits + 2), 0ULL);
    fib[1] = 1ULL;
    for (int i = 2; i <= bits + 1; ++i) {
        fib[static_cast<std::size_t>(i)] =
            fib[static_cast<std::size_t>(i - 1)] + fib[static_cast<std::size_t>(i - 2)];
    }

    // Count of n in [1, 2^exponent] with no adjacent 1 bits.
    return fib[static_cast<std::size_t>(bits + 1)];
}

u64 brute_small(const int exponent) {
    const u64 limit = 1ULL << exponent;
    u64 count = 0;
    for (u64 n = 1; n <= limit; ++n) {
        if ((n ^ (2ULL * n) ^ (3ULL * n)) == 0ULL) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(3) != brute_small(3)) {
        std::cerr << "Checkpoint failed for exponent=3" << '\n';
        return false;
    }
    if (solve(10) != brute_small(10)) {
        std::cerr << "Checkpoint failed for exponent=10" << '\n';
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
