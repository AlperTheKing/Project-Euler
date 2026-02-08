#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 10'000'000;
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
    std::vector<bool> is_composite(static_cast<std::size_t>(limit + 1), false);
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (!is_composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            if (i <= limit / i) {
                for (int j = i * i; j <= limit; j += i) {
                    is_composite[static_cast<std::size_t>(j)] = true;
                }
            }
        }
    }
    return primes;
}

u64 max_for_pair(const int p, const int q, const int limit) {
    if (static_cast<u64>(p) * static_cast<u64>(q) > static_cast<u64>(limit)) {
        return 0ULL;
    }
    u64 best = 0ULL;
    for (u64 pa = static_cast<u64>(p); pa <= static_cast<u64>(limit / q);) {
        u64 value = pa * static_cast<u64>(q);
        while (value <= static_cast<u64>(limit)) {
            best = std::max(best, value);
            if (value > static_cast<u64>(limit / q)) {
                break;
            }
            value *= static_cast<u64>(q);
        }
        if (pa > static_cast<u64>(limit / p)) {
            break;
        }
        pa *= static_cast<u64>(p);
    }
    return best;
}

u64 solve(const int limit) {
    const std::vector<int> primes = sieve_primes(limit);
    std::vector<unsigned char> seen(static_cast<std::size_t>(limit + 1), 0U);

    u64 sum = 0ULL;
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const int p = primes[i];
        if (static_cast<u64>(p) * 2ULL > static_cast<u64>(limit)) {
            break;
        }
        for (std::size_t j = i + 1; j < primes.size(); ++j) {
            const int q = primes[j];
            if (static_cast<u64>(p) * static_cast<u64>(q) > static_cast<u64>(limit)) {
                break;
            }
            const u64 m = max_for_pair(p, q, limit);
            if (m > 0ULL && seen[static_cast<std::size_t>(m)] == 0U) {
                seen[static_cast<std::size_t>(m)] = 1U;
                sum += m;
            }
        }
    }
    return sum;
}

u64 brute_solve(const int limit) {
    std::unordered_map<u64, int> best_by_pair;
    for (int n = 2; n <= limit; ++n) {
        int x = n;
        std::vector<int> distinct_primes;
        for (int p = 2; p <= x / p; ++p) {
            if (x % p != 0) {
                continue;
            }
            distinct_primes.push_back(p);
            while (x % p == 0) {
                x /= p;
            }
            if (distinct_primes.size() > 2U) {
                break;
            }
        }
        if (distinct_primes.size() <= 2U && x > 1) {
            distinct_primes.push_back(x);
        }
        if (distinct_primes.size() == 2U) {
            int p = distinct_primes[0];
            int q = distinct_primes[1];
            if (p > q) {
                std::swap(p, q);
            }
            const u64 key = (static_cast<u64>(static_cast<std::uint32_t>(p)) << 32U) |
                            static_cast<u64>(static_cast<std::uint32_t>(q));
            auto it = best_by_pair.find(key);
            if (it == best_by_pair.end() || n > it->second) {
                best_by_pair[key] = n;
            }
        }
    }

    std::vector<unsigned char> seen(static_cast<std::size_t>(limit + 1), 0U);
    u64 sum = 0ULL;
    for (const auto& entry : best_by_pair) {
        const int value = entry.second;
        if (seen[static_cast<std::size_t>(value)] == 0U) {
            seen[static_cast<std::size_t>(value)] = 1U;
            sum += static_cast<u64>(value);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (max_for_pair(2, 3, 100) != 96ULL) {
        std::cerr << "Checkpoint failed for M(2,3,100)" << '\n';
        return false;
    }
    if (max_for_pair(3, 5, 100) != 75ULL) {
        std::cerr << "Checkpoint failed for M(3,5,100)" << '\n';
        return false;
    }
    if (max_for_pair(2, 73, 100) != 0ULL) {
        std::cerr << "Checkpoint failed for M(2,73,100)" << '\n';
        return false;
    }
    if (solve(100) != 2262ULL) {
        std::cerr << "Checkpoint failed for S(100)" << '\n';
        return false;
    }
    if (solve(300) != brute_solve(300)) {
        std::cerr << "Checkpoint failed for brute-force cross-check at N=300" << '\n';
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
