#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Options {
    int limit = 100000;
    int index = 10000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--index=", options.index)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 1 && options.index >= 1 && options.index <= options.limit;
}

int solve(const int limit, const int index) {
    std::vector<int> rad(static_cast<std::size_t>(limit + 1), 1);
    std::vector<char> is_prime(static_cast<std::size_t>(limit + 1), true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }

    for (int p = 2; p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int m = p; m <= limit; m += p) {
            rad[static_cast<std::size_t>(m)] *= p;
            if (m > p) {
                is_prime[static_cast<std::size_t>(m)] = false;
            }
        }
    }

    std::vector<std::pair<int, int>> values;
    values.reserve(static_cast<std::size_t>(limit));
    for (int n = 1; n <= limit; ++n) {
        values.push_back({rad[static_cast<std::size_t>(n)], n});
    }

    std::sort(values.begin(), values.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) {
            return a.first < b.first;
        }
        return a.second < b.second;
    });

    return values[static_cast<std::size_t>(index - 1)].second;
}

bool run_checkpoints() {
    if (solve(10, 4) != 8) {
        std::cerr << "Checkpoint failed for limit=10,index=4" << '\n';
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

    std::cout << solve(options.limit, options.index) << '\n';
    return 0;
}
