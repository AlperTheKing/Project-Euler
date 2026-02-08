#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    int n = 1000;
    bool run_checkpoints = true;
};

constexpr u64 kMod = 10000000000ULL;

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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

u64 mod_pow(u64 base, int exponent) {
    u64 result = 1ULL;
    base %= kMod;

    while (exponent > 0) {
        if ((exponent & 1) != 0) {
            result = static_cast<u64>((__uint128_t)result * base % kMod);
        }
        base = static_cast<u64>((__uint128_t)base * base % kMod);
        exponent >>= 1;
    }

    return result;
}

u64 solve(const int n) {
    u64 total = 0ULL;
    for (int k = 1; k <= n; ++k) {
        total = (total + mod_pow(static_cast<u64>(k), k)) % kMod;
    }
    return total;
}

bool run_checkpoints() {
    if (solve(10) != 405071317ULL) {
        std::cerr << "Checkpoint failed for n=10" << '\n';
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

    std::cout << solve(options.n) << '\n';
    return 0;
}
