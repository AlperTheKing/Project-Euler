#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int limit = 100000000;
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
    return options.limit >= 4;
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    if (limit >= 0) {
        is_prime[0] = 0U;
    }
    if (limit >= 1) {
        is_prime[1] = 0U;
    }
    for (int p = 2; static_cast<i64>(p) * p <= limit; ++p) {
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

i64 solve(const int limit) {
    const int max_q = (limit - 1) / 2;
    const std::vector<int> primes = sieve_primes(max_q);

    i64 count = 0;
    const int total = static_cast<int>(primes.size());
    for (int i = 0; i < total; ++i) {
        const i64 p = primes[static_cast<std::size_t>(i)];
        if (p * p >= limit) {
            break;
        }
        const i64 q_max = (limit - 1) / p;
        const auto it = std::upper_bound(primes.begin() + i, primes.end(), static_cast<int>(q_max));
        count += static_cast<i64>(it - (primes.begin() + i));
    }
    return count;
}

i64 brute_small(const int limit) {
    int count = 0;
    for (int n = 4; n < limit; ++n) {
        int x = n;
        int omega = 0;
        for (int p = 2; p * p <= x; ++p) {
            while (x % p == 0) {
                x /= p;
                ++omega;
                if (omega > 2) {
                    break;
                }
            }
            if (omega > 2) {
                break;
            }
        }
        if (x > 1) {
            ++omega;
        }
        if (omega == 2) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(30) != 10) {
        std::cerr << "Checkpoint failed for limit 30 sample" << '\n';
        return false;
    }
    if (solve(10000) != brute_small(10000)) {
        std::cerr << "Checkpoint failed for brute cross-check at limit 10000" << '\n';
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
