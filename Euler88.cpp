#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int k_max = 12000;
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
        if (parse_int_after_prefix(arg, "--k-max=", options.k_max)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.k_max >= 2;
}

void dfs(const int start,
         const i64 product,
         const i64 sum,
         const int factors,
         const int k_max,
         std::vector<i64>& best) {
    const i64 k = product - sum + factors;
    if (k > k_max) {
        return;
    }

    if (factors >= 2) {
        if (product < best[static_cast<std::size_t>(k)]) {
            best[static_cast<std::size_t>(k)] = product;
        }
    }

    for (int x = start; ; ++x) {
        const i64 next_product = product * x;
        if (next_product > 2LL * k_max) {
            break;
        }
        dfs(x, next_product, sum + x, factors + 1, k_max, best);
    }
}

i64 solve(const int k_max) {
    std::vector<i64> best(static_cast<std::size_t>(k_max + 1), std::numeric_limits<i64>::max());
    dfs(2, 1, 0, 0, k_max, best);

    std::set<i64> unique_values;
    for (int k = 2; k <= k_max; ++k) {
        unique_values.insert(best[static_cast<std::size_t>(k)]);
    }

    i64 total = 0;
    for (const i64 x : unique_values) {
        total += x;
    }
    return total;
}

bool run_checkpoints() {
    if (solve(6) != 30LL) {
        std::cerr << "Checkpoint failed for k_max=6" << '\n';
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

    std::cout << solve(options.k_max) << '\n';
    return 0;
}
