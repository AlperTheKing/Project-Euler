#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = long long;

struct Options {
    u64 limit = 10000000000000000ULL;
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

std::vector<int> primes_below_100() {
    constexpr int max_n = 100;
    std::array<bool, max_n + 1> composite{};
    std::vector<int> primes;
    for (int i = 2; i < max_n; ++i) {
        if (!composite[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
        for (int p : primes) {
            const int v = i * p;
            if (v >= max_n) {
                break;
            }
            composite[static_cast<std::size_t>(v)] = true;
            if (i % p == 0) {
                break;
            }
        }
    }
    return primes;
}

i64 binom(const int n, const int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    i64 numer = 1;
    i64 denom = 1;
    for (int i = 1; i <= k; ++i) {
        numer *= (n - k + i);
        denom *= i;
    }
    return numer / denom;
}

u64 solve(const u64 limit) {
    if (limit <= 1) {
        return 0;
    }

    const std::vector<int> primes = primes_below_100();
    std::array<u64, 26> subset_sums{};

    const u64 max_value = limit - 1;

    const auto dfs = [&](const auto& self, const std::size_t index, const int taken, const u64 product) -> void {
        if (taken >= 4) {
            subset_sums[static_cast<std::size_t>(taken)] += max_value / product;
        }

        if (index >= primes.size()) {
            return;
        }

        for (std::size_t i = index; i < primes.size(); ++i) {
            const u64 p = static_cast<u64>(primes[i]);
            if (product > max_value / p) {
                continue;
            }
            self(self, i + 1, taken + 1, product * p);
        }
    };

    dfs(dfs, 0U, 0, 1ULL);

    i64 answer = 0;
    for (int k = 4; k <= static_cast<int>(primes.size()); ++k) {
        const i64 coeff = ((k - 4) % 2 == 0 ? 1 : -1) * binom(k - 1, 3);
        answer += coeff * static_cast<i64>(subset_sums[static_cast<std::size_t>(k)]);
    }

    return static_cast<u64>(answer);
}

u64 brute(const int limit) {
    const std::vector<int> primes = primes_below_100();
    u64 count = 0;
    for (int x = 1; x < limit; ++x) {
        int distinct = 0;
        int t = x;
        for (int p : primes) {
            if (t % p == 0) {
                ++distinct;
                while (t % p == 0) {
                    t /= p;
                }
            }
        }
        if (distinct >= 4) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(1000ULL) != 23ULL) {
        std::cerr << "Checkpoint failed for limit 1000" << '\n';
        return false;
    }
    if (solve(100000ULL) != brute(100000)) {
        std::cerr << "Checkpoint failed for brute-force cross-check at 100000" << '\n';
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
