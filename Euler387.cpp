#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 limit = 100'000'000'000'000ULL;
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 100ULL;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

bool is_prime(const u64 n) {
    if (n < 2ULL) {
        return false;
    }
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0ULL) {
            return false;
        }
    }

    u64 d = n - 1ULL;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    for (u64 a : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL}) {
        if (a >= n) {
            continue;
        }
        u64 x = pow_mod(a, d, n);
        if (x == 1ULL || x == n - 1ULL) {
            continue;
        }
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1ULL) {
                witness = false;
                break;
            }
        }
        if (witness) {
            return false;
        }
    }
    return true;
}

u64 dfs_harshad(const u64 n, const int digit_sum, const u64 harshad_limit, const u64 prime_limit) {
    u64 acc = 0ULL;

    if (n >= 10ULL && n % static_cast<u64>(digit_sum) == 0ULL) {
        const u64 q = n / static_cast<u64>(digit_sum);
        if (is_prime(q)) {
            for (int d : {1, 3, 7, 9}) {
                const u64 cand = n * 10ULL + static_cast<u64>(d);
                if (cand < prime_limit && is_prime(cand)) {
                    acc += cand;
                }
            }
        }
    }

    if (n > harshad_limit / 10ULL) {
        return acc;
    }

    for (int d = 0; d <= 9; ++d) {
        const u64 next = n * 10ULL + static_cast<u64>(d);
        const int next_sum = digit_sum + d;
        if (next_sum != 0 && next % static_cast<u64>(next_sum) == 0ULL) {
            acc += dfs_harshad(next, next_sum, harshad_limit, prime_limit);
        }
    }
    return acc;
}

u64 solve(const u64 limit) {
    const u64 harshad_limit = (limit - 1ULL) / 10ULL;
    u64 total = 0ULL;
    for (int d = 1; d <= 9; ++d) {
        total += dfs_harshad(static_cast<u64>(d), d, harshad_limit, limit);
    }
    return total;
}

u64 brute_small(const u64 limit) {
    auto digit_sum = [](u64 n) {
        int s = 0;
        while (n > 0ULL) {
            s += static_cast<int>(n % 10ULL);
            n /= 10ULL;
        }
        return s;
    };
    auto is_right_trunc_harshad = [&](u64 n) {
        while (n > 0ULL) {
            const int s = digit_sum(n);
            if (s == 0 || n % static_cast<u64>(s) != 0ULL) {
                return false;
            }
            n /= 10ULL;
        }
        return true;
    };

    u64 sum = 0ULL;
    for (u64 p = 2ULL; p < limit; ++p) {
        if (!is_prime(p)) {
            continue;
        }
        const u64 h = p / 10ULL;
        if (h < 10ULL) {
            continue;
        }
        if (!is_right_trunc_harshad(h)) {
            continue;
        }
        const int s = digit_sum(h);
        if (h % static_cast<u64>(s) == 0ULL && is_prime(h / static_cast<u64>(s))) {
            sum += p;
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(10'000ULL) != 90'619ULL) {
        std::cerr << "Checkpoint failed: limit=10000 sample" << '\n';
        return false;
    }
    if (solve(1'000'000ULL) != brute_small(1'000'000ULL)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for limit=1e6" << '\n';
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
    std::cout << solve(options.limit) << '\n';
    return 0;
}
