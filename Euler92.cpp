#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 10000000;
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

int next_square_digit_sum(int n) {
    int sum = 0;
    while (n > 0) {
        const int d = n % 10;
        sum += d * d;
        n /= 10;
    }
    return sum;
}

int chain_end(int n, std::vector<int>& memo) {
    std::vector<int> path;
    while (n >= static_cast<int>(memo.size()) || memo[static_cast<std::size_t>(n)] == 0) {
        path.push_back(n);
        n = next_square_digit_sum(n);
    }

    const int end = memo[static_cast<std::size_t>(n)];
    for (int x : path) {
        if (x < static_cast<int>(memo.size())) {
            memo[static_cast<std::size_t>(x)] = end;
        }
    }
    return end;
}

int solve(const int limit) {
    const int max_sum = 9 * 9 * 7;  // for numbers below 10^7
    std::vector<int> memo(static_cast<std::size_t>(max_sum + 1), 0);
    memo[1] = 1;
    memo[89] = 89;

    int count = 0;
    for (int n = 1; n < limit; ++n) {
        const int s = next_square_digit_sum(n);
        if (chain_end(s, memo) == 89) {
            ++count;
        }
    }

    return count;
}

bool run_checkpoints() {
    if (solve(100) != 80) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
        return false;
    }
    if (solve(10) != 7) {
        std::cerr << "Checkpoint failed for limit=10" << '\n';
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
