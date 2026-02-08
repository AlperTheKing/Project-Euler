#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int index = 30;
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
        if (parse_int_after_prefix(arg, "--index=", options.index)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.index >= 1;
}

int digit_sum(u64 n) {
    int sum = 0;
    while (n > 0) {
        sum += static_cast<int>(n % 10ULL);
        n /= 10ULL;
    }
    return sum;
}

int digits_count(u64 n) {
    int count = 0;
    while (n > 0) {
        ++count;
        n /= 10ULL;
    }
    return std::max(1, count);
}

std::vector<u64> generate_terms(const u64 limit) {
    const int max_digit_sum = 9 * digits_count(limit);
    std::set<u64> values;

    for (int s = 2; s <= max_digit_sum; ++s) {
        u128 value = static_cast<u128>(s) * static_cast<u128>(s);
        while (value <= limit) {
            const u64 v = static_cast<u64>(value);
            if (v >= 10ULL && digit_sum(v) == s) {
                values.insert(v);
            }
            value *= static_cast<u128>(s);
        }
    }

    return std::vector<u64>(values.begin(), values.end());
}

u64 solve(const int index) {
    u64 limit = 1000000ULL;

    for (;;) {
        std::vector<u64> terms = generate_terms(limit);
        if (static_cast<int>(terms.size()) >= index) {
            return terms[static_cast<std::size_t>(index - 1)];
        }
        limit *= 10ULL;
    }
}

bool run_checkpoints() {
    if (solve(2) != 512ULL) {
        std::cerr << "Checkpoint failed for second term" << '\n';
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

    std::cout << solve(options.index) << '\n';
    return 0;
}
