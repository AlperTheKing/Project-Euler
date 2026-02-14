#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

namespace {

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
    for (const char c : tail) {
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

std::vector<bool> prime_sieve(const int limit) {
    std::vector<bool> prime(static_cast<std::size_t>(limit + 1), true);
    if (limit >= 0) {
        prime[0] = false;
    }
    if (limit >= 1) {
        prime[1] = false;
    }
    for (int p = 2; p <= limit / p; ++p) {
        if (!prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            prime[static_cast<std::size_t>(q)] = false;
        }
    }
    return prime;
}

bool is_prime_trial(const int n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (int p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

bool is_circular_prime(const int n, const std::vector<bool>& sieve) {
    std::string s = std::to_string(n);

    for (std::size_t shift = 0; shift < s.size(); ++shift) {
        std::rotate(s.begin(), s.begin() + 1, s.end());
        const int rot = std::stoi(s);

        bool prime = false;
        if (rot < static_cast<int>(sieve.size())) {
            prime = sieve[static_cast<std::size_t>(rot)];
        } else {
            prime = is_prime_trial(rot);
        }
        if (!prime) {
            return false;
        }
    }

    return true;
}

int solve(const int limit) {
    const auto sieve = prime_sieve(limit + 100);

    int count = 0;
    for (int n = 2; n < limit; ++n) {
        if (!sieve[static_cast<std::size_t>(n)]) {
            continue;
        }
        if (is_circular_prime(n, sieve)) {
            ++count;
        }
    }

    return count;
}

bool run_checkpoints() {
    const auto sieve = prime_sieve(1000);
    if (!is_circular_prime(197, sieve)) {
        std::cerr << "Checkpoint failed for 197" << '\n';
        return false;
    }
    if (solve(100) != 13) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
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
