#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 n = 1'000'000'000'000'000ULL;
    u64 mod = 1'000'000'000ULL;
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
        if (parse_u64_after_prefix(arg, "--n=", options.n) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL && options.mod >= 2ULL;
}

u64 sum_squares_mod(u64 n, const u64 mod) {
    if (n == 0ULL) {
        return 0ULL;
    }
    u64 a = n;
    u64 b = n + 1ULL;
    u64 c = 2ULL * n + 1ULL;

    if ((a & 1ULL) == 0ULL) {
        a /= 2ULL;
    } else {
        b /= 2ULL;
    }
    if (a % 3ULL == 0ULL) {
        a /= 3ULL;
    } else if (b % 3ULL == 0ULL) {
        b /= 3ULL;
    } else {
        c /= 3ULL;
    }

    const __uint128_t x = static_cast<__uint128_t>(a % mod) * static_cast<__uint128_t>(b % mod);
    const __uint128_t y = (x % mod) * static_cast<__uint128_t>(c % mod);
    return static_cast<u64>(y % mod);
}

u64 sum_squares_range_mod(const u64 l, const u64 r, const u64 mod) {
    const u64 right = sum_squares_mod(r, mod);
    const u64 left = sum_squares_mod(l - 1ULL, mod);
    return (right + mod - left) % mod;
}

u64 solve(const u64 n, const u64 mod) {
    u64 ans = 0ULL;
    for (u64 l = 1ULL; l <= n;) {
        const u64 q = n / l;
        const u64 r = n / q;
        const u64 sum_d2 = sum_squares_range_mod(l, r, mod);
        const __uint128_t term = static_cast<__uint128_t>(q % mod) * static_cast<__uint128_t>(sum_d2);
        ans = (ans + static_cast<u64>(term % mod)) % mod;
        l = r + 1ULL;
    }
    return ans;
}

u64 brute_sigma2_sum(const u64 n) {
    u64 s = 0ULL;
    for (u64 x = 1ULL; x <= n; ++x) {
        u64 sigma2 = 0ULL;
        for (u64 d = 1ULL; d * d <= x; ++d) {
            if (x % d != 0ULL) {
                continue;
            }
            sigma2 += d * d;
            const u64 e = x / d;
            if (e != d) {
                sigma2 += e * e;
            }
        }
        s += sigma2;
    }
    return s;
}

bool run_checkpoints() {
    const u64 mod = 1'000'000'000ULL;
    if (solve(1ULL, mod) != 1ULL ||
        solve(2ULL, mod) != 6ULL ||
        solve(3ULL, mod) != 16ULL ||
        solve(4ULL, mod) != 37ULL ||
        solve(5ULL, mod) != 63ULL ||
        solve(6ULL, mod) != 113ULL) {
        std::cerr << "Checkpoint failed: first SIGMA2 values" << '\n';
        return false;
    }
    if (solve(1000ULL, mod) != brute_sigma2_sum(1000ULL) % mod) {
        std::cerr << "Checkpoint failed: brute-force cross-check for n=1000" << '\n';
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
