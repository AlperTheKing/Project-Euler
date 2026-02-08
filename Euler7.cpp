#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int index = 10001;
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

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--index=", options.index)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.index > 0;
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit) + 1ULL, true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }

    for (int p = 2; p * p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            primes.push_back(p);
        }
    }

    return primes;
}

int solve(const int index) {
    if (index == 1) {
        return 2;
    }

    int upper = 64;
    if (index >= 6) {
        const double n = static_cast<double>(index);
        upper = static_cast<int>(n * (std::log(n) + std::log(std::log(n)))) + 16;
    }

    while (true) {
        const std::vector<int> primes = sieve_primes(upper);
        if (static_cast<int>(primes.size()) >= index) {
            return primes[static_cast<std::size_t>(index - 1)];
        }
        upper *= 2;
    }
}

bool run_checkpoints() {
    if (solve(6) != 13) {
        std::cerr << "Checkpoint failed for index=6" << '\n';
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

    std::cout << solve(options.index) << '\n';
    return 0;
}
