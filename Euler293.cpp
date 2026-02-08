#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 1000000000ULL;
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
    return options.limit >= 3ULL;
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        exp >>= 1U;
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
        d >>= 1U;
        ++s;
    }

    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0ULL) {
            continue;
        }
        u64 x = mod_pow(a, d, n);
        if (x == 1ULL || x == n - 1ULL) {
            continue;
        }
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = static_cast<u64>((static_cast<u128>(x) * x) % n);
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

std::vector<int> prime_prefix() {
    return {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
}

void generate_admissible(const std::vector<int>& primes,
                         const int idx,
                         const u64 current,
                         const u64 limit,
                         std::set<u64>& out) {
    if (idx >= static_cast<int>(primes.size())) {
        return;
    }
    const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(idx)]);
    u64 value = current;
    while (value <= (limit - 1ULL) / p) {
        value *= p;
        out.insert(value);
        generate_admissible(primes, idx + 1, value, limit, out);
    }
}

int pseudo_fortunate(const u64 n) {
    for (int m = 3;; m += 2) {
        if (is_prime(n + static_cast<u64>(m))) {
            return m;
        }
    }
}

u64 solve(const u64 limit) {
    std::set<u64> admissible;
    const auto primes = prime_prefix();
    generate_admissible(primes, 0, 1ULL, limit, admissible);

    std::set<int> pseudo_values;
    for (u64 n : admissible) {
        pseudo_values.insert(pseudo_fortunate(n));
    }

    u64 sum = 0ULL;
    for (int m : pseudo_values) {
        sum += static_cast<u64>(m);
    }
    return sum;
}

bool is_admissible_bruteforce(u64 n) {
    if ((n & 1ULL) == 1ULL) {
        return false;
    }
    std::vector<u64> factors;
    u64 x = n;
    for (u64 p = 2ULL; p * p <= x; ++p) {
        if (x % p != 0ULL) {
            continue;
        }
        factors.push_back(p);
        while (x % p == 0ULL) {
            x /= p;
        }
    }
    if (x > 1ULL) {
        factors.push_back(x);
    }
    if (factors.empty() || factors[0] != 2ULL) {
        return false;
    }
    u64 candidate = 2ULL;
    std::size_t idx = 0;
    while (idx < factors.size()) {
        if (factors[idx] != candidate) {
            return false;
        }
        ++idx;
        ++candidate;
        while (!is_prime(candidate)) {
            ++candidate;
        }
    }
    return true;
}

u64 solve_bruteforce(const u64 limit) {
    std::set<int> pseudo_values;
    for (u64 n = 2; n < limit; n += 2) {
        if (!is_admissible_bruteforce(n)) {
            continue;
        }
        pseudo_values.insert(pseudo_fortunate(n));
    }
    u64 sum = 0ULL;
    for (int m : pseudo_values) {
        sum += static_cast<u64>(m);
    }
    return sum;
}

bool run_checkpoints() {
    if (pseudo_fortunate(16ULL) != 3) {
        std::cerr << "Checkpoint failed for N=16 sample" << '\n';
        return false;
    }
    if (pseudo_fortunate(630ULL) != 11) {
        std::cerr << "Checkpoint failed for N=630 sample" << '\n';
        return false;
    }
    if (solve(10000ULL) != solve_bruteforce(10000ULL)) {
        std::cerr << "Checkpoint failed for brute cross-check at limit=10000" << '\n';
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
