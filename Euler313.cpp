#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int prime_limit = 1000000;
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
        if (parse_int_after_prefix(arg, "--prime-limit=", options.prime_limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.prime_limit >= 2;
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;
    for (int p = 2; static_cast<u64>(p) * p <= static_cast<u64>(limit); ++p) {
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

u64 count_grids_for_prime_square(const u64 p2) {
    if ((p2 & 1ULL) == 0ULL) {
        return 0ULL;
    }
    const u64 h = (p2 + 13ULL) / 2ULL;
    const u64 low = h / 4ULL + 1ULL;
    const u64 high = (h >= 2ULL) ? ((h - 2ULL) / 3ULL) : 0ULL;
    if (high < low) {
        return 0ULL;
    }
    return 2ULL * (high - low + 1ULL);
}

u64 solve(const int prime_limit) {
    const std::vector<int> primes = sieve_primes(prime_limit - 1);
    u64 total = 0ULL;
    for (int p : primes) {
        const u64 p2 = static_cast<u64>(p) * static_cast<u64>(p);
        total += count_grids_for_prime_square(p2);
    }
    return total;
}

bool run_checkpoints() {
    if (solve(100) != 5482ULL) {
        std::cerr << "Checkpoint failed for prime limit 100 sample" << '\n';
        return false;
    }
    if (solve(10) != 8ULL) {
        std::cerr << "Checkpoint failed for prime limit 10" << '\n';
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
    std::cout << solve(options.prime_limit) << '\n';
    return 0;
}
