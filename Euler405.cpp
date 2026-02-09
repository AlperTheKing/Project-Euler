#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
constexpr u64 kMod = 410338673ULL;   // 17^7
constexpr u64 kPhi = 386201104ULL;   // phi(17^7)

struct Options {
    i64 k = 1000000000000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoll(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_i64_after_prefix(arg, "--k=", options.k)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.k >= 0;
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0) {
        if (e & 1ULL) {
            result = static_cast<u64>((__uint128_t)result * cur % mod);
        }
        cur = static_cast<u64>((__uint128_t)cur * cur % mod);
        e >>= 1ULL;
    }
    return result;
}

u64 f_mod_n(const u64 n, const u64 mod) {
    // Closed form:
    // f(n) = 1 - (-1)^n / 15 - (4/3) * 2^n + (2/5) * 4^n.
    const u64 inv3 = mod_pow(3, mod - mod / 17 - 1, mod);
    const u64 inv5 = mod_pow(5, mod - mod / 17 - 1, mod);
    const u64 inv15 = static_cast<u64>((__uint128_t)inv3 * inv5 % mod);

    const u64 sign = (n & 1ULL) ? (mod - 1ULL) : 1ULL;
    const u64 p2 = mod_pow(2, n, mod);
    const u64 p4 = mod_pow(4, n, mod);

    u64 ans = 1ULL;
    ans = (ans + mod - static_cast<u64>((__uint128_t)sign * inv15 % mod)) % mod;
    ans = (ans + mod - static_cast<u64>((__uint128_t)4ULL * p2 % mod * inv3 % mod)) % mod;
    ans = (ans + static_cast<u64>((__uint128_t)2ULL * p4 % mod * inv5 % mod)) % mod;
    return ans;
}

u64 solve(i64 k) {
    // n = 10^k, and we need f(n) mod 17^7.
    // For 2^n and 4^n modulo 17^7, reduce exponent modulo phi(17^7).
    const u64 n_mod_phi = mod_pow(10, static_cast<u64>(k), kPhi);

    // For k>=1, n is even => (-1)^n = 1. For k=0, n=1 is odd.
    if (k == 0) {
        return f_mod_n(1ULL, kMod);
    }

    // Use n_mod_phi for exponentiations; parity is known even.
    const u64 inv3 = mod_pow(3, kPhi - 1, kMod);
    const u64 inv5 = mod_pow(5, kPhi - 1, kMod);
    const u64 inv15 = static_cast<u64>((__uint128_t)inv3 * inv5 % kMod);
    const u64 p2 = mod_pow(2, n_mod_phi, kMod);
    const u64 p4 = mod_pow(4, n_mod_phi, kMod);

    u64 ans = 1ULL;
    ans = (ans + kMod - inv15) % kMod;  // (-1)^n = 1
    ans = (ans + kMod - static_cast<u64>((__uint128_t)4ULL * p2 % kMod * inv3 % kMod)) % kMod;
    ans = (ans + static_cast<u64>((__uint128_t)2ULL * p4 % kMod * inv5 % kMod)) % kMod;
    return ans;
}

bool run_checkpoints() {
    if (f_mod_n(1ULL, kMod) != 0ULL) {
        std::cerr << "Checkpoint failed: f(1)\n";
        return false;
    }
    if (f_mod_n(4ULL, kMod) != 82ULL) {
        std::cerr << "Checkpoint failed: f(4)\n";
        return false;
    }
    if (f_mod_n(1000000000ULL, kMod) != 126897180ULL) {
        std::cerr << "Checkpoint failed: f(10^9) mod 17^7\n";
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

    std::cout << solve(options.k) << '\n';
    return 0;
}
