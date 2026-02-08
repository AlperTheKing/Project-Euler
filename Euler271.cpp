#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
    u64 n = 13082761331670030ULL;
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n > 2;
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<i128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<i128>(cur) * cur) % mod);
        exp >>= 1U;
    }
    return result;
}

std::vector<u64> distinct_prime_factors(u64 n) {
    std::vector<u64> factors;
    if ((n & 1ULL) == 0ULL) {
        factors.push_back(2ULL);
        while ((n & 1ULL) == 0ULL) {
            n >>= 1U;
        }
    }
    for (u64 p = 3; p * p <= n; p += 2ULL) {
        if (n % p != 0ULL) {
            continue;
        }
        factors.push_back(p);
        while (n % p == 0ULL) {
            n /= p;
        }
    }
    if (n > 1ULL) {
        factors.push_back(n);
    }
    return factors;
}

i64 extended_gcd(const i64 a, const i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

u64 mod_inverse(const u64 a, const u64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    if (g != 1) {
        return 0;
    }
    i64 r = x % static_cast<i64>(mod);
    if (r < 0) {
        r += static_cast<i64>(mod);
    }
    return static_cast<u64>(r);
}

u64 crt_pair(const u64 a1, const u64 m1, const u64 a2, const u64 m2) {
    const u64 inv = mod_inverse(m1 % m2, m2);
    const u64 t = static_cast<u64>((static_cast<i128>((a2 + m2 - (a1 % m2)) % m2) * inv) % m2);
    return static_cast<u64>(a1 + static_cast<i128>(m1) * t);
}

u64 solve(const u64 n) {
    const std::vector<u64> primes = distinct_prime_factors(n);

    std::vector<std::vector<u64>> roots_by_prime;
    roots_by_prime.reserve(primes.size());
    for (u64 p : primes) {
        std::vector<u64> roots;
        for (u64 x = 1; x < p; ++x) {
            if (mod_pow(x, 3, p) == 1ULL) {
                roots.push_back(x);
            }
        }
        roots_by_prime.push_back(std::move(roots));
    }

    u64 sum = 0;
    const auto dfs = [&](auto&& self, std::size_t idx, u64 residue, u64 modulus) -> void {
        if (idx == primes.size()) {
            if (residue > 1ULL && residue < n) {
                sum += residue;
            }
            return;
        }
        const u64 p = primes[idx];
        for (u64 root : roots_by_prime[idx]) {
            const u64 next_residue = crt_pair(residue, modulus, root, p);
            self(self, idx + 1, next_residue, modulus * p);
        }
    };
    dfs(dfs, 0, 0ULL, 1ULL);
    return sum;
}

bool run_checkpoints() {
    if (solve(91) != 363ULL) {
        std::cerr << "Checkpoint failed for n=91 sample" << '\n';
        return false;
    }
    if (solve(7) != 6ULL) {
        std::cerr << "Checkpoint failed for n=7" << '\n';
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
