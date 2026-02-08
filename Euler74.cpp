#include <array>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Options {
    int limit = 1000000;
    int target_length = 60;
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
        if (parse_int_after_prefix(arg, "--target-length=", options.target_length)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 1 && options.target_length >= 1;
}

int next_value(int n, const std::array<int, 10>& fac) {
    int sum = 0;
    do {
        sum += fac[static_cast<std::size_t>(n % 10)];
        n /= 10;
    } while (n > 0);
    return sum;
}

int chain_length(int start,
                 const std::array<int, 10>& fac,
                 std::unordered_map<int, int>& memo) {
    std::vector<int> order;
    order.reserve(80U);
    std::unordered_map<int, int> seen;

    int n = start;
    while (true) {
        auto it_memo = memo.find(n);
        if (it_memo != memo.end()) {
            int len = it_memo->second;
            for (int i = static_cast<int>(order.size()) - 1; i >= 0; --i) {
                ++len;
                memo[order[static_cast<std::size_t>(i)]] = len;
            }
            return memo[start];
        }

        auto it_seen = seen.find(n);
        if (it_seen != seen.end()) {
            const int loop_start = it_seen->second;
            const int loop_len = static_cast<int>(order.size()) - loop_start;

            for (int i = loop_start; i < static_cast<int>(order.size()); ++i) {
                memo[order[static_cast<std::size_t>(i)]] = loop_len;
            }

            int len = loop_len;
            for (int i = loop_start - 1; i >= 0; --i) {
                ++len;
                memo[order[static_cast<std::size_t>(i)]] = len;
            }
            return memo[start];
        }

        seen[n] = static_cast<int>(order.size());
        order.push_back(n);
        n = next_value(n, fac);
    }
}

int solve(const int limit, const int target_length) {
    std::array<int, 10> fac{};
    fac[0] = 1;
    for (int d = 1; d <= 9; ++d) {
        fac[static_cast<std::size_t>(d)] = fac[static_cast<std::size_t>(d - 1)] * d;
    }

    std::unordered_map<int, int> memo;
    memo.reserve(2000000U);

    int count = 0;
    for (int n = 1; n < limit; ++n) {
        if (chain_length(n, fac, memo) == target_length) {
            ++count;
        }
    }

    return count;
}

bool run_checkpoints() {
    std::array<int, 10> fac{};
    fac[0] = 1;
    for (int d = 1; d <= 9; ++d) {
        fac[static_cast<std::size_t>(d)] = fac[static_cast<std::size_t>(d - 1)] * d;
    }

    std::unordered_map<int, int> memo;
    if (chain_length(69, fac, memo) != 5) {
        std::cerr << "Checkpoint failed for start=69" << '\n';
        return false;
    }
    if (chain_length(78, fac, memo) != 4) {
        std::cerr << "Checkpoint failed for start=78" << '\n';
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

    std::cout << solve(options.limit, options.target_length) << '\n';
    return 0;
}
