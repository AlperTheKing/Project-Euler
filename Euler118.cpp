#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using i64 = std::int64_t;

struct Options {
    int max_digit = 9;
    bool run_checkpoints = true;
};

struct PrimeEntry {
    u32 value = 0;
    int mask = 0;
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
        if (parse_int_after_prefix(arg, "--max-digit=", options.max_digit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_digit >= 1 && options.max_digit <= 9;
}

bool is_prime(u32 n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2U) == 0U) {
        return n == 2;
    }
    for (u32 p = 3; static_cast<std::uint64_t>(p) * p <= n; p += 2) {
        if ((n % p) == 0U) {
            return false;
        }
    }
    return true;
}

void generate_primes_rec(const int max_digit,
                         int used_mask,
                         u32 value,
                         std::vector<PrimeEntry>& entries) {
    if (value > 0 && is_prime(value)) {
        entries.push_back({value, used_mask});
    }

    for (int d = 1; d <= max_digit; ++d) {
        const int bit = 1 << (d - 1);
        if ((used_mask & bit) != 0) {
            continue;
        }
        generate_primes_rec(max_digit, used_mask | bit, value * 10U + static_cast<u32>(d), entries);
    }
}

i64 count_sets(const int max_digit) {
    std::vector<PrimeEntry> entries;
    entries.reserve(50000);
    generate_primes_rec(max_digit, 0, 0, entries);

    std::sort(entries.begin(), entries.end(), [](const PrimeEntry& a, const PrimeEntry& b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }
        return a.mask < b.mask;
    });
    entries.erase(std::unique(entries.begin(), entries.end(), [](const PrimeEntry& a, const PrimeEntry& b) {
                      return a.value == b.value && a.mask == b.mask;
                  }),
                  entries.end());

    const int full_mask = (1 << max_digit) - 1;
    std::vector<i64> memo(static_cast<std::size_t>(1 << max_digit), -1);

    const auto dfs = [&](auto&& self, int remaining_mask) -> i64 {
        if (remaining_mask == 0) {
            return 1;
        }
        i64& memo_cell = memo[static_cast<std::size_t>(remaining_mask)];
        if (memo_cell != -1) {
            return memo_cell;
        }

        const int anchor_bit = remaining_mask & -remaining_mask;
        i64 total = 0;

        for (const PrimeEntry& e : entries) {
            if ((e.mask & anchor_bit) == 0) {
                continue;
            }
            if ((e.mask & remaining_mask) != e.mask) {
                continue;
            }
            total += self(self, remaining_mask ^ e.mask);
        }

        memo_cell = total;
        return total;
    };

    return dfs(dfs, full_mask);
}

bool run_checkpoints() {
    if (count_sets(4) != 9) {
        std::cerr << "Checkpoint failed for digits 1..4" << '\n';
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

    std::cout << count_sets(options.max_digit) << '\n';
    return 0;
}
