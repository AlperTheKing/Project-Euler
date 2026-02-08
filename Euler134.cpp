#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 1000000;
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
    return options.limit >= 5;
}

i64 extended_gcd(i64 a, i64 b, i64& x, i64& y) {
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

i64 mod_inverse(i64 a, i64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(a, mod, x, y);
    if (g != 1) {
        return -1;
    }
    x %= mod;
    if (x < 0) {
        x += mod;
    }
    return x;
}

i64 pair_value(const int p1, const int p2) {
    i64 mod = 1;
    int x = p1;
    while (x > 0) {
        mod *= 10;
        x /= 10;
    }

    const i64 inv = mod_inverse(p2 % mod, mod);
    const i64 k = (static_cast<i64>(p1) * inv) % mod;
    return static_cast<i64>(p2) * k;
}

std::vector<int> primes_up_to(const int limit) {
    std::vector<bool> sieve(static_cast<std::size_t>(limit + 1), true);
    sieve[0] = false;
    sieve[1] = false;

    for (int i = 2; static_cast<i64>(i) * i <= limit; ++i) {
        if (!sieve[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            sieve[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (sieve[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

i64 solve(const int limit) {
    const std::vector<int> primes = primes_up_to(limit + 200000);

    i64 sum = 0;
    for (std::size_t i = 0; i + 1 < primes.size(); ++i) {
        const int p1 = primes[i];
        const int p2 = primes[i + 1];
        if (p1 < 5) {
            continue;
        }
        if (p1 >= limit) {
            break;
        }
        sum += pair_value(p1, p2);
    }

    return sum;
}

bool run_checkpoints() {
    if (pair_value(19, 23) != 1219) {
        std::cerr << "Checkpoint failed for pair (19,23)" << '\n';
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
