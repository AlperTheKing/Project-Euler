#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int limit = 100000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 2;
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, u64 mod) {
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

int multiplicative_order_10(const int p) {
    int order = p - 1;
    int x = order;

    std::vector<int> factors;
    for (int f = 2; static_cast<std::int64_t>(f) * f <= x; ++f) {
        if (x % f != 0) {
            continue;
        }
        factors.push_back(f);
        while (x % f == 0) {
            x /= f;
        }
    }
    if (x > 1) {
        factors.push_back(x);
    }

    for (int q : factors) {
        while ((order % q) == 0 && pow_mod(10ULL, static_cast<u64>(order / q), static_cast<u64>(p)) == 1ULL) {
            order /= q;
        }
    }

    return order;
}

std::vector<int> primes_below(const int limit) {
    std::vector<bool> sieve(static_cast<std::size_t>(limit), true);
    if (limit > 0) {
        sieve[0] = false;
    }
    if (limit > 1) {
        sieve[1] = false;
    }

    for (int i = 2; static_cast<std::int64_t>(i) * i < limit; ++i) {
        if (!sieve[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j < limit; j += i) {
            sieve[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i < limit; ++i) {
        if (sieve[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

bool can_divide_some_R10n(const int p) {
    if (p == 2 || p == 3 || p == 5) {
        return false;
    }

    int ord = multiplicative_order_10(p);
    while ((ord % 2) == 0) {
        ord /= 2;
    }
    while ((ord % 5) == 0) {
        ord /= 5;
    }

    return ord == 1;
}

std::int64_t solve(const int limit) {
    const std::vector<int> primes = primes_below(limit);

    std::int64_t sum = 0;
    for (int p : primes) {
        if (!can_divide_some_R10n(p)) {
            sum += p;
        }
    }

    return sum;
}

bool run_checkpoints() {
    if (solve(100) != 918) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
        return false;
    }
    if (!can_divide_some_R10n(17) || can_divide_some_R10n(19) || can_divide_some_R10n(3)) {
        std::cerr << "Checkpoint failed for prime behavior (3, 17, 19)" << '\n';
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
