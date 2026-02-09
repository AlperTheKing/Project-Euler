#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int n = 1'000'000;
    u64 mod = 1'000'000'007ULL;
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
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1 && options.mod >= 2ULL;
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

u64 mod_inverse(const u64 a, const u64 mod) {
    return mod_pow(a % mod, mod - 2ULL, mod);
}

u64 solve(const int n, const u64 mod) {
    u64 ans = 0ULL;
    for (int k = 1; k <= n; ++k) {
        const i64 kk = static_cast<i64>(k);
        i64 t = 1LL - kk * kk;
        t %= static_cast<i64>(mod);
        if (t < 0) {
            t += static_cast<i64>(mod);
        }
        const u64 tm = static_cast<u64>(t);

        u64 contribution = 0ULL;
        if (tm == 1ULL) {
            contribution = static_cast<u64>(n) % mod;
        } else {
            const u64 pn = mod_pow(tm, static_cast<u64>(n), mod);
            const u64 numerator = static_cast<u64>((static_cast<u128>(tm) * ((pn + mod - 1ULL) % mod)) % mod);
            const u64 denom = (tm + mod - 1ULL) % mod;
            contribution = static_cast<u64>((static_cast<u128>(numerator) * mod_inverse(denom, mod)) % mod);
        }
        ans += contribution;
        if (ans >= mod) {
            ans -= mod;
        }
    }
    return ans;
}

u64 brute(const int n, const u64 mod) {
    u64 s = 0ULL;
    for (int k = 1; k <= n; ++k) {
        const i64 kk = static_cast<i64>(k);
        i64 u = 1LL - kk * kk;
        u %= static_cast<i64>(mod);
        if (u < 0) {
            u += static_cast<i64>(mod);
        }
        u64 term = 1ULL;
        for (int p = 1; p <= n; ++p) {
            term = static_cast<u64>((static_cast<u128>(term) * static_cast<u64>(u)) % mod);
            s += term;
            if (s >= mod) {
                s -= mod;
            }
        }
    }
    return s;
}

bool run_checkpoints() {
    const u64 mod = 1'000'000'007ULL;
    if (solve(4, mod) != 51'160ULL) {
        std::cerr << "Checkpoint failed: S(4)=51160" << '\n';
        return false;
    }
    if (solve(20, mod) != brute(20, mod)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for n=20" << '\n';
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
    std::cout << solve(options.n, options.mod) << '\n';
    return 0;
}
