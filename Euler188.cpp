#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 base = 1777;
    int height = 1855;
    u64 mod = 100000000ULL;
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
        if (parse_u64_after_prefix(arg, "--base=", options.base) ||
            parse_int_after_prefix(arg, "--height=", options.height) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.height >= 1 && options.mod >= 1;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    u64 result = 1 % mod;
    base %= mod;

    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }

    return result;
}

u64 phi(u64 n) {
    u64 result = n;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        while (n % p == 0) {
            n /= p;
        }
        result = result / p * (p - 1);
    }
    if (n > 1) {
        result = result / n * (n - 1);
    }
    return result;
}

u64 tetration_mod(const u64 base, const int height, const u64 mod) {
    if (mod == 1) {
        return 0;
    }
    if (height == 1) {
        return base % mod;
    }

    const u64 ph = phi(mod);
    u64 exp = tetration_mod(base, height - 1, ph);
    if (std::gcd(base, mod) != 1ULL) {
        exp += ph;
    }

    return pow_mod(base, exp, mod);
}

bool run_checkpoints() {
    if (tetration_mod(3, 3, 100) != 87ULL) {
        std::cerr << "Checkpoint failed for 3^^3 mod 100" << '\n';
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

    std::cout << tetration_mod(options.base, options.height, options.mod) << '\n';
    return 0;
}
