#include <cstdint>
#include <iostream>
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

    return options.limit >= 3;
}

std::vector<bool> sieve(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    is_prime[0] = false;
    is_prime[1] = false;

    for (int p = 2; p <= limit / p; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    return is_prime;
}

int solve(const int limit) {
    const std::vector<bool> is_prime = sieve(limit - 1);

    std::vector<int> primes;
    for (int n = 2; n < limit; ++n) {
        if (is_prime[static_cast<std::size_t>(n)]) {
            primes.push_back(n);
        }
    }

    std::vector<i64> prefix(primes.size() + 1U, 0LL);
    for (std::size_t i = 0; i < primes.size(); ++i) {
        prefix[i + 1U] = prefix[i] + primes[i];
    }

    int best_prime = 0;
    int best_len = 0;

    for (std::size_t i = 0; i < primes.size(); ++i) {
        for (std::size_t j = i + static_cast<std::size_t>(best_len) + 1U; j <= primes.size(); ++j) {
            const i64 sum = prefix[j] - prefix[i];
            if (sum >= limit) {
                break;
            }
            if (is_prime[static_cast<std::size_t>(sum)]) {
                best_len = static_cast<int>(j - i);
                best_prime = static_cast<int>(sum);
            }
        }
    }

    return best_prime;
}

bool run_checkpoints() {
    if (solve(100) != 41) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
        return false;
    }
    if (solve(1000) != 953) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
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
