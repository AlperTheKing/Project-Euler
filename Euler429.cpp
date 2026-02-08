#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int n = 100'000'000;
    u64 mod = 1'000'000'009ULL;
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

std::vector<int> primes_up_to(const int n) {
    std::vector<bool> is_composite(static_cast<std::size_t>(n + 1), false);
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            if (i <= n / i) {
                for (int j = i * i; j <= n; j += i) {
                    is_composite[static_cast<std::size_t>(j)] = true;
                }
            }
        }
    }
    return primes;
}

u64 exponent_in_factorial(const int n, const int p) {
    u64 e = 0ULL;
    int x = n;
    while (x > 0) {
        x /= p;
        e += static_cast<u64>(x);
    }
    return e;
}

u64 solve(const int n, const u64 mod) {
    const std::vector<int> primes = primes_up_to(n);
    u64 ans = 1ULL;
    for (int p : primes) {
        const u64 e = exponent_in_factorial(n, p);
        const u64 pe2 = mod_pow(static_cast<u64>(p), 2ULL * e, mod);
        const u64 factor = (1ULL + pe2) % mod;
        ans = static_cast<u64>((static_cast<u128>(ans) * factor) % mod);
    }
    return ans;
}

u64 brute_unitary_sum_square_factorial(const int n) {
    u64 fact = 1ULL;
    for (int i = 2; i <= n; ++i) {
        fact *= static_cast<u64>(i);
    }
    u64 sum = 0ULL;
    for (u64 d = 1ULL; d <= fact; ++d) {
        if (fact % d != 0ULL) {
            continue;
        }
        const u64 q = fact / d;
        u64 a = d;
        u64 b = q;
        while (b != 0ULL) {
            const u64 t = a % b;
            a = b;
            b = t;
        }
        if (a == 1ULL) {
            sum += d * d;
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(4, 1'000'000'009ULL) != 650ULL) {
        std::cerr << "Checkpoint failed: S(4!)=650" << '\n';
        return false;
    }
    const u64 mod = 1'000'000'009ULL;
    if (solve(8, mod) != brute_unitary_sum_square_factorial(8) % mod) {
        std::cerr << "Checkpoint failed: brute-force cross-check for n=8" << '\n';
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
