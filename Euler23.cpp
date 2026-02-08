#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 28123;
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

    i64 parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
        if (parsed > static_cast<i64>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    value = static_cast<int>(parsed);
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
    return options.limit >= 1;
}

std::vector<int> proper_divisor_sums(const int limit) {
    std::vector<int> sums(static_cast<std::size_t>(limit + 1), 0);
    for (int d = 1; d <= limit / 2; ++d) {
        for (int m = d * 2; m <= limit; m += d) {
            sums[static_cast<std::size_t>(m)] += d;
        }
    }
    return sums;
}

i64 solve(const int limit) {
    const std::vector<int> sums = proper_divisor_sums(limit);

    std::vector<int> abundant;
    abundant.reserve(static_cast<std::size_t>(limit / 2));
    for (int n = 12; n <= limit; ++n) {
        if (sums[static_cast<std::size_t>(n)] > n) {
            abundant.push_back(n);
        }
    }

    std::vector<bool> is_abundant_sum(static_cast<std::size_t>(limit + 1), false);
    for (std::size_t i = 0; i < abundant.size(); ++i) {
        for (std::size_t j = i; j < abundant.size(); ++j) {
            const int value = abundant[i] + abundant[j];
            if (value > limit) {
                break;
            }
            is_abundant_sum[static_cast<std::size_t>(value)] = true;
        }
    }

    i64 total = 0;
    for (int n = 1; n <= limit; ++n) {
        if (!is_abundant_sum[static_cast<std::size_t>(n)]) {
            total += n;
        }
    }

    return total;
}

bool run_checkpoints() {
    if (solve(50) != 891) {
        std::cerr << "Checkpoint failed for limit=50" << '\n';
        return false;
    }
    if (solve(100) != 2766) {
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
