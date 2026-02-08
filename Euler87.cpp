#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 50000000;
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

std::vector<int> sieve_primes(const int limit) {
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
    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            primes.push_back(p);
        }
    }
    return primes;
}

int solve(const int limit) {
    const int p_limit = static_cast<int>(std::sqrt(limit)) + 1;
    const std::vector<int> primes = sieve_primes(p_limit);

    std::vector<i64> square;
    std::vector<i64> cube;
    std::vector<i64> fourth;

    for (int p : primes) {
        const i64 p2 = static_cast<i64>(p) * p;
        if (p2 < limit) {
            square.push_back(p2);
        }

        const i64 p3 = p2 * p;
        if (p3 < limit) {
            cube.push_back(p3);
        }

        const i64 p4 = p3 * p;
        if (p4 < limit) {
            fourth.push_back(p4);
        }
    }

    std::unordered_set<int> values;
    values.reserve(2000000U);

    for (i64 a : square) {
        for (i64 b : cube) {
            if (a + b >= limit) {
                break;
            }
            for (i64 c : fourth) {
                const i64 s = a + b + c;
                if (s >= limit) {
                    break;
                }
                values.insert(static_cast<int>(s));
            }
        }
    }

    return static_cast<int>(values.size());
}

bool run_checkpoints() {
    if (solve(50) != 4) {
        std::cerr << "Checkpoint failed for limit=50" << '\n';
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
