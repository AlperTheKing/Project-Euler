#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int limit = 1'000'000;
    u64 mod = 1'000'000'000ULL;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
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
    u64 parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2 && options.mod >= 2ULL;
}

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return result;
}

u64 f_value(const u64 n, const u64 mod) {
    u64 x = 0ULL;
    std::vector<u64> seen;
    seen.reserve(32);

    for (int iter = 0; iter < 80; ++iter) {
        const u64 nx = mod_pow(n, x, mod);
        if (nx == x) {
            return (x > 0ULL) ? x : 0ULL;
        }
        bool repeated = false;
        for (const u64 v : seen) {
            if (v == nx) {
                repeated = true;
                break;
            }
        }
        if (repeated) {
            return 0ULL;
        }
        seen.push_back(nx);
        x = nx;
    }
    const u64 nx = mod_pow(n, x, mod);
    if (nx == x && x > 0ULL) {
        return x;
    }
    return 0ULL;
}

u64 solve(const int limit, const u64 mod) {
    u64 sum = 0ULL;
    for (int n = 2; n <= limit; ++n) {
        sum += f_value(static_cast<u64>(n), mod);
    }
    return sum;
}

bool run_checkpoints() {
    const u64 mod = 1'000'000'000ULL;
    if (f_value(4ULL, mod) != 411'728'896ULL) {
        std::cerr << "Checkpoint failed: f(4)" << '\n';
        return false;
    }
    if (f_value(10ULL, mod) != 0ULL) {
        std::cerr << "Checkpoint failed: f(10)=0" << '\n';
        return false;
    }
    if (f_value(157ULL, mod) != 743'757ULL) {
        std::cerr << "Checkpoint failed: f(157)" << '\n';
        return false;
    }
    if (solve(1'000, mod) != 442'530'011'399ULL) {
        std::cerr << "Checkpoint failed: sample sum for n<=1000" << '\n';
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
    std::cout << solve(options.limit, options.mod) << '\n';
    return 0;
}
