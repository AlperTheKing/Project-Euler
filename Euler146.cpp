#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 150000000ULL;
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 10;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    base %= mod;
    u64 result = 1 % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

bool is_prime(const u64 n) {
    if (n < 2) {
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

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    static constexpr std::array<u64, 7> kBases{
        2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};

    for (u64 a : kBases) {
        if (a % n == 0ULL) {
            continue;
        }
        u64 x = pow_mod(a, d, n);
        if (x == 1ULL || x == n - 1) {
            continue;
        }
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
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

u64 solve(const u64 limit) {
    static constexpr std::array<u64, 6> kGood{1ULL, 3ULL, 7ULL, 9ULL, 13ULL, 27ULL};
    static constexpr std::array<u64, 8> kBad{5ULL, 11ULL, 15ULL, 17ULL, 19ULL, 21ULL, 23ULL, 25ULL};

    u64 sum = 0;
    for (u64 n = 10; n < limit; n += 10) {
        const u64 square = n * n;
        if (square % 3ULL == 0ULL || square % 7ULL == 0ULL || square % 13ULL == 0ULL) {
            continue;
        }

        bool ok = true;
        for (u64 k : kGood) {
            if (!is_prime(square + k)) {
                ok = false;
                break;
            }
        }
        if (!ok) {
            continue;
        }

        for (u64 k : kBad) {
            if (is_prime(square + k)) {
                ok = false;
                break;
            }
        }
        if (ok) {
            sum += n;
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (!is_prime(2ULL) || !is_prime(3ULL) || is_prime(1ULL) || is_prime(21ULL)) {
        std::cerr << "Checkpoint failed for primality tester" << '\n';
        return false;
    }
    if (solve(1000000ULL) != 1242490ULL) {
        std::cerr << "Checkpoint failed for limit 1,000,000" << '\n';
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
