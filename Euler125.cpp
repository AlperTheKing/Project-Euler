#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    i64 limit = 100000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2;
}

bool is_palindrome(i64 n) {
    std::string s = std::to_string(n);
    std::string t = s;
    std::reverse(t.begin(), t.end());
    return s == t;
}

i64 solve(const i64 limit) {
    int max_base = 1;
    while (static_cast<i64>(max_base) * max_base < limit) {
        ++max_base;
    }

    std::vector<i64> prefix(static_cast<std::size_t>(max_base + 1), 0);
    for (int i = 1; i <= max_base; ++i) {
        prefix[static_cast<std::size_t>(i)] = prefix[static_cast<std::size_t>(i - 1)] +
                                              static_cast<i64>(i) * static_cast<i64>(i);
    }

    std::set<i64> values;

    for (int start = 1; start <= max_base; ++start) {
        for (int end = start + 1; end <= max_base; ++end) {
            const i64 sum = prefix[static_cast<std::size_t>(end)] - prefix[static_cast<std::size_t>(start - 1)];
            if (sum >= limit) {
                break;
            }
            if (is_palindrome(sum)) {
                values.insert(sum);
            }
        }
    }

    i64 total = 0;
    for (i64 v : values) {
        total += v;
    }
    return total;
}

bool run_checkpoints() {
    if (solve(1000) != 4164) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
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
