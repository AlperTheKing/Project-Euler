#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 start = 100000000000000ULL;
    int count = 100000;
    u64 modulo = 1234567891011ULL;
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
        if (parse_u64_after_prefix(arg, "--start=", options.start) ||
            parse_int_after_prefix(arg, "--count=", options.count) ||
            parse_u64_after_prefix(arg, "--modulo=", options.modulo)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.start >= 2ULL && options.count >= 1 && options.modulo >= 2ULL;
}

std::vector<int> sieve_base_primes(const int limit) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;
    for (int p = 2; static_cast<u64>(p) * static_cast<u64>(p) <= static_cast<u64>(limit); ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = 0U;
        }
    }
    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] != 0U) {
            primes.push_back(p);
        }
    }
    return primes;
}

std::vector<u64> first_primes_after(u64 start, const int count) {
    const int base_limit = 20000000;  // sqrt(1e14 + margin) is around 1e7.
    const std::vector<int> base_primes = sieve_base_primes(base_limit);

    const u64 segment_size = 4000000ULL;
    std::vector<u64> out;
    out.reserve(static_cast<std::size_t>(count));

    u64 low = start + 1ULL;
    while (static_cast<int>(out.size()) < count) {
        const u64 high = low + segment_size - 1ULL;
        std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(segment_size), 1U);

        for (int p : base_primes) {
            const u64 pp = static_cast<u64>(p);
            if (pp * pp > high) {
                break;
            }
            u64 first = (low + pp - 1ULL) / pp * pp;
            if (first < pp * pp) {
                first = pp * pp;
            }
            for (u64 x = first; x <= high; x += pp) {
                is_prime[static_cast<std::size_t>(x - low)] = 0U;
            }
        }

        for (u64 i = 0; i < segment_size && static_cast<int>(out.size()) < count; ++i) {
            const u64 value = low + i;
            if (value >= 2ULL && is_prime[static_cast<std::size_t>(i)] != 0U) {
                out.push_back(value);
            }
        }
        low = high + 1ULL;
    }
    return out;
}

std::pair<u64, u64> fib_pair_mod(const u64 n, const u64 mod) {
    if (n == 0ULL) {
        return {0ULL, 1ULL % mod};
    }
    const auto [a, b] = fib_pair_mod(n >> 1U, mod);
    const u64 two_b = static_cast<u64>((2ULL * b) % mod);
    const u64 c = static_cast<u64>((static_cast<u128>(a) * ((two_b + mod - a) % mod)) % mod);
    const u64 d = static_cast<u64>((static_cast<u128>(a) * a + static_cast<u128>(b) * b) % mod);
    if ((n & 1ULL) == 0ULL) {
        return {c, d};
    }
    return {d, static_cast<u64>((c + d) % mod)};
}

u64 solve(const u64 start, const int count, const u64 mod) {
    const std::vector<u64> primes = first_primes_after(start, count);
    u64 sum = 0ULL;
    for (u64 p : primes) {
        sum += fib_pair_mod(p, mod).first;
        if (sum >= mod) {
            sum %= mod;
        }
    }
    return sum % mod;
}

u64 solve_bruteforce(const u64 start, const int count, const u64 mod) {
    auto is_prime_trial = [](u64 x) {
        if (x < 2ULL) {
            return false;
        }
        if ((x & 1ULL) == 0ULL) {
            return x == 2ULL;
        }
        for (u64 d = 3ULL; d * d <= x; d += 2ULL) {
            if (x % d == 0ULL) {
                return false;
            }
        }
        return true;
    };

    std::vector<u64> primes;
    for (u64 x = start + 1ULL; static_cast<int>(primes.size()) < count; ++x) {
        if (is_prime_trial(x)) {
            primes.push_back(x);
        }
    }
    u64 sum = 0ULL;
    for (u64 p : primes) {
        sum = (sum + fib_pair_mod(p, mod).first) % mod;
    }
    return sum;
}

bool run_checkpoints() {
    if (fib_pair_mod(10ULL, 1000000007ULL).first != 55ULL) {
        std::cerr << "Checkpoint failed for Fibonacci doubling F(10)" << '\n';
        return false;
    }
    if (solve(1000ULL, 30, 1000000007ULL) != solve_bruteforce(1000ULL, 30, 1000000007ULL)) {
        std::cerr << "Checkpoint failed for brute cross-check on small prime range" << '\n';
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
    std::cout << solve(options.start, options.count, options.modulo) << '\n';
    return 0;
}
