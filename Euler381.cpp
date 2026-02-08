#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 100'000'000;
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
    return options.limit >= 6;
}

std::vector<int> primes_below(const int limit) {
    std::vector<bool> is_composite(static_cast<std::size_t>(limit), false);
    std::vector<int> primes;
    for (int i = 2; i < limit; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            if (i <= (limit - 1) / i) {
                for (int j = i * i; j < limit; j += i) {
                    is_composite[static_cast<std::size_t>(j)] = true;
                }
            }
        }
    }
    return primes;
}

u64 inv8_mod_prime(const u64 p) {
    const u64 r = p & 7ULL;
    const u64 k = 8ULL - r;
    return (k * p + 1ULL) / 8ULL;
}

u64 S_of_prime(const u64 p) {
    const u64 inv8 = inv8_mod_prime(p);
    return (p + p - ((3ULL * inv8) % p)) % p;
}

u64 brute_S_of_prime(const u64 p) {
    u64 fact = 1ULL;
    for (u64 i = 1; i < p; ++i) {
        fact = (fact * i) % p;
    }
    u64 term = fact;
    u64 sum = 0ULL;
    for (u64 k = 1; k <= 5; ++k) {
        sum = (sum + term) % p;
        const u64 denom = p - k;
        u64 inv = 1ULL;
        for (u64 x = 1; x < p; ++x) {
            if ((denom * x) % p == 1ULL) {
                inv = x;
                break;
            }
        }
        term = (term * inv) % p;
    }
    return sum;
}

u64 solve(const int limit) {
    const std::vector<int> primes = primes_below(limit);
    u64 sum = 0ULL;
    for (int p : primes) {
        if (p < 5) {
            continue;
        }
        sum += S_of_prime(static_cast<u64>(p));
    }
    return sum;
}

bool run_checkpoints() {
    if (S_of_prime(7ULL) != 4ULL) {
        std::cerr << "Checkpoint failed: S(7)" << '\n';
        return false;
    }
    for (u64 p : {5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL}) {
        if (S_of_prime(p) != brute_S_of_prime(p)) {
            std::cerr << "Checkpoint failed: direct factorial check for prime " << p << '\n';
            return false;
        }
    }
    if (solve(100) != 480ULL) {
        std::cerr << "Checkpoint failed: sum S(p), 5<=p<100" << '\n';
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
