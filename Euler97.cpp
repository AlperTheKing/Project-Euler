#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 coefficient = 28433ULL;
    u64 base = 2ULL;
    u64 exponent = 7830457ULL;
    u64 increment = 1ULL;
    u64 mod = 10000000000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--coefficient=", options.coefficient) ||
            parse_u64_after_prefix(arg, "--base=", options.base) ||
            parse_u64_after_prefix(arg, "--exponent=", options.exponent) ||
            parse_u64_after_prefix(arg, "--increment=", options.increment) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.mod > 0;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exponent, const u64 mod) {
    u64 result = 1ULL % mod;
    base %= mod;

    while (exponent > 0) {
        if ((exponent & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exponent >>= 1ULL;
    }

    return result;
}

u64 solve(const u64 coefficient, const u64 base, const u64 exponent, const u64 increment, const u64 mod) {
    const u64 power = pow_mod(base, exponent, mod);
    const u64 scaled = mul_mod(coefficient % mod, power, mod);
    return (scaled + increment % mod) % mod;
}

bool run_checkpoints() {
    if (pow_mod(2ULL, 10ULL, 1000ULL) != 24ULL) {
        std::cerr << "Checkpoint failed for modular exponentiation" << '\n';
        return false;
    }

    const u64 small = solve(28433ULL, 2ULL, 20ULL, 1ULL, 10000000000ULL);
    if (small != 9814161409ULL) {
        std::cerr << "Checkpoint failed for reduced exponent scenario" << '\n';
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

    std::cout << solve(options.coefficient, options.base, options.exponent, options.increment, options.mod) << '\n';
    return 0;
}
