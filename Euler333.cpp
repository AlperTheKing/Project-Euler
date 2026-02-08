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
    return options.prime_limit >= 3;
}

std::vector<std::vector<int>> build_terms(const int limit) {
    std::vector<std::vector<int>> rows;
    for (u64 p2 = 1ULL; p2 <= static_cast<u64>(limit); p2 <<= 1ULL) {
        std::vector<int> row;
        for (u64 value = p2; value <= static_cast<u64>(limit); value *= 3ULL) {
            row.push_back(static_cast<int>(value));
            if (value > static_cast<u64>(limit) / 3ULL) {
                break;
            }
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

std::vector<std::uint32_t> partition_counts(const int limit) {
    const std::vector<std::vector<int>> rows = build_terms(limit);
    std::vector<std::uint32_t> count(static_cast<std::size_t>(limit + 1), 0U);

    const auto dfs = [&](auto&& self, const int i, const int j_limit, const int sum) -> void {
        if (i == static_cast<int>(rows.size())) {
            if (sum > 0) {
                ++count[static_cast<std::size_t>(sum)];
            }
            return;
        }

        self(self, i + 1, j_limit, sum);  // Skip this 2^i layer.

        const int upper = std::min(j_limit, static_cast<int>(rows[static_cast<std::size_t>(i)].size()));
        for (int j = 0; j < upper; ++j) {
            const int value = rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            if (sum + value > limit) {
                break;
            }
            self(self, i + 1, j, sum + value);
        }
    };
    dfs(dfs, 0, 1000, 0);
    return count;
}

u64 solve(const int prime_limit) {
    const int limit = prime_limit - 1;
    const std::vector<std::uint32_t> counts = partition_counts(limit);

    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(prime_limit), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;
    for (int p = 2; static_cast<long long>(p) * p < prime_limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (int q = p * p; q < prime_limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = 0U;
        }
    }

    u64 sum = 0ULL;
    for (int q = 2; q < prime_limit; ++q) {
        if (is_prime[static_cast<std::size_t>(q)] != 0U &&
            counts[static_cast<std::size_t>(q)] == 1U) {
            sum += static_cast<u64>(q);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(100) != 233ULL) {
        std::cerr << "Checkpoint failed for stated sample q<100" << '\n';
        return false;
    }
    const auto counts = partition_counts(100);
    if (counts[11] != 2U || counts[17] != 1U) {
        std::cerr << "Checkpoint failed for P(11)=2 and P(17)=1 examples" << '\n';
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
