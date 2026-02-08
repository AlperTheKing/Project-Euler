#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 250250;
    int divisor = 250;
    u64 modulo = 10000000000000000ULL;
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
        parsed = parsed * 10 + static_cast<u64>(c - '0');
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--divisor=", options.divisor) ||
            parse_u64_after_prefix(arg, "--mod=", options.modulo)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 1 && options.divisor >= 1 && options.modulo > 0;
}

int powmod_int(int base, int exp, const int mod) {
    int result = 1 % mod;
    int cur = base % mod;

    while (exp > 0) {
        if (exp & 1) {
            result = static_cast<int>((static_cast<long long>(result) * cur) % mod);
        }
        cur = static_cast<int>((static_cast<long long>(cur) * cur) % mod);
        exp >>= 1;
    }

    return result;
}

u64 solve(const int limit, const int divisor, const u64 modulo) {
    std::vector<u64> dp(static_cast<std::size_t>(divisor), 0);
    std::vector<u64> next(static_cast<std::size_t>(divisor), 0);
    dp[0] = 1;

    for (int i = 1; i <= limit; ++i) {
        const int residue = powmod_int(i, i, divisor);

        for (int r = 0; r < divisor; ++r) {
            next[static_cast<std::size_t>(r)] = dp[static_cast<std::size_t>(r)];
        }

        for (int r = 0; r < divisor; ++r) {
            const int to = (r + residue) % divisor;
            u64& cell = next[static_cast<std::size_t>(to)];
            cell += dp[static_cast<std::size_t>(r)];
            if (cell >= modulo) {
                cell -= modulo;
            }
        }

        dp.swap(next);
    }

    const u64 all_divisible = dp[0];
    return (all_divisible + modulo - 1) % modulo;  // Exclude empty subset.
}

bool run_checkpoints() {
    if (solve(10, 250, 10000000000000000ULL) != 5ULL) {
        std::cerr << "Checkpoint failed for first 10 values" << '\n';
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

    std::cout << solve(options.limit, options.divisor, options.modulo) << '\n';
    return 0;
}
