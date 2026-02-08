#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int prime_limit = 5000;
    u64 modulo = 10000000000000000ULL;
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
        if (parse_int_after_prefix(arg, "--prime-limit=", options.prime_limit) ||
            parse_u64_after_prefix(arg, "--mod=", options.modulo)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.prime_limit >= 2 && options.modulo > 0;
}

std::vector<int> primes_below(const int limit) {
    std::vector<std::uint8_t> sieve(static_cast<std::size_t>(limit), 1);
    sieve[0] = 0;
    sieve[1] = 0;

    for (int i = 2; i * i < limit; ++i) {
        if (!sieve[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j < limit; j += i) {
            sieve[static_cast<std::size_t>(j)] = 0;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i < limit; ++i) {
        if (sieve[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

u64 solve(const int prime_limit, const u64 modulo) {
    const std::vector<int> primes = primes_below(prime_limit);

    int max_sum = 0;
    for (int p : primes) {
        max_sum += p;
    }

    std::vector<u64> ways(static_cast<std::size_t>(max_sum + 1), 0);
    ways[0] = 1;

    int current_sum = 0;
    for (int p : primes) {
        for (int s = current_sum; s >= 0; --s) {
            const u64 add = ways[static_cast<std::size_t>(s)];
            if (add == 0) {
                continue;
            }
            u64& dst = ways[static_cast<std::size_t>(s + p)];
            dst += add;
            if (dst >= modulo) {
                dst -= modulo;
            }
        }
        current_sum += p;
    }

    std::vector<std::uint8_t> is_prime_sum(static_cast<std::size_t>(max_sum + 1), 1);
    is_prime_sum[0] = 0;
    is_prime_sum[1] = 0;
    for (int i = 2; i * i <= max_sum; ++i) {
        if (!is_prime_sum[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= max_sum; j += i) {
            is_prime_sum[static_cast<std::size_t>(j)] = 0;
        }
    }

    u64 answer = 0;
    for (int s = 2; s <= max_sum; ++s) {
        if (!is_prime_sum[static_cast<std::size_t>(s)]) {
            continue;
        }
        answer += ways[static_cast<std::size_t>(s)];
        answer %= modulo;
    }

    return answer;
}

bool run_checkpoints() {
    if (solve(10, 10000000000000000ULL) != 7ULL) {
        std::cerr << "Checkpoint failed for primes below 10" << '\n';
        return false;
    }
    if (solve(20, 10000000000000000ULL) != 65ULL) {
        std::cerr << "Checkpoint failed for primes below 20" << '\n';
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

    std::cout << solve(options.prime_limit, options.modulo) << '\n';
    return 0;
}
