#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int target = 200;
    u64 initial_limit = 1000000000ULL;
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
        parsed = parsed * 10 + static_cast<u64>(c - '0');
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
        if (parse_int_after_prefix(arg, "--target=", options.target) ||
            parse_u64_after_prefix(arg, "--initial-limit=", options.initial_limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.target >= 1 && options.initial_limit >= 100;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    u64 result = 1 % mod;
    base %= mod;

    while (exp > 0) {
        if (exp & 1ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1;
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
        if (n % p == 0) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1;
        ++s;
    }

    constexpr u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (a % n == 0) {
            continue;
        }
        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) {
            continue;
        }

        bool witnessed = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                witnessed = false;
                break;
            }
        }
        if (witnessed) {
            return false;
        }
    }

    return true;
}

std::vector<int> primes_up_to(const int n) {
    std::vector<std::uint8_t> is_prime_small(static_cast<std::size_t>(n + 1), 1);
    is_prime_small[0] = 0;
    is_prime_small[1] = 0;

    for (int i = 2; i * i <= n; ++i) {
        if (!is_prime_small[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime_small[static_cast<std::size_t>(j)] = 0;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime_small[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

bool contains_200(u64 x) {
    while (x >= 200) {
        if (x % 1000ULL == 200ULL) {
            return true;
        }
        x /= 10ULL;
    }
    return false;
}

bool is_prime_proof(const u64 x) {
    u64 pow10 = 1;
    while (pow10 <= x) {
        const int digit = static_cast<int>((x / pow10) % 10ULL);
        const u64 base = x - static_cast<u64>(digit) * pow10;

        for (int repl = 0; repl <= 9; ++repl) {
            if (repl == digit) {
                continue;
            }
            if (pow10 * 10ULL > x && repl == 0) {
                continue;
            }

            const u64 candidate = base + static_cast<u64>(repl) * pow10;
            if (is_prime(candidate)) {
                return false;
            }
        }

        pow10 *= 10ULL;
    }

    return true;
}

std::vector<u64> generate_squbes(const u64 limit) {
    const u64 max_p = static_cast<u64>(std::sqrt(static_cast<long double>(limit / 8ULL))) + 5ULL;
    const std::vector<int> primes = primes_up_to(static_cast<int>(max_p));

    std::vector<u64> squbes;
    for (int p : primes) {
        const u128 p2 = static_cast<u128>(p) * static_cast<u128>(p);
        if (p2 * 8U > limit) {
            break;
        }

        for (int q : primes) {
            if (q == p) {
                continue;
            }

            const u128 q3 = static_cast<u128>(q) * q * q;
            const u128 value = p2 * q3;
            if (value > limit) {
                break;
            }

            squbes.push_back(static_cast<u64>(value));
        }
    }

    std::sort(squbes.begin(), squbes.end());
    squbes.erase(std::unique(squbes.begin(), squbes.end()), squbes.end());
    return squbes;
}

std::vector<u64> first_squbes(const std::size_t k) {
    u64 limit = 1000;
    while (true) {
        std::vector<u64> squbes = generate_squbes(limit);
        if (squbes.size() >= k) {
            squbes.resize(k);
            return squbes;
        }
        limit *= 2;
    }
}

u64 solve(const int target, u64 limit) {
    while (true) {
        const std::vector<u64> squbes = generate_squbes(limit);

        std::vector<u64> matches;
        matches.reserve(squbes.size() / 4);

        for (u64 x : squbes) {
            if (!contains_200(x)) {
                continue;
            }
            if (is_prime_proof(x)) {
                matches.push_back(x);
            }
        }

        if (static_cast<int>(matches.size()) >= target) {
            return matches[static_cast<std::size_t>(target - 1)];
        }

        limit *= 2ULL;
    }
}

bool run_checkpoints() {
    const std::vector<u64> first = first_squbes(5);
    const std::vector<u64> expected = {72ULL, 108ULL, 200ULL, 392ULL, 500ULL};
    if (first != expected) {
        std::cerr << "Checkpoint failed for first squbes" << '\n';
        return false;
    }

    if (solve(2, 2000000ULL) != 1992008ULL) {
        std::cerr << "Checkpoint failed for second prime-proof sqube containing 200" << '\n';
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

    std::cout << solve(options.target, options.initial_limit) << '\n';
    return 0;
}
