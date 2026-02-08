#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int max_index = 5000;
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
        if (parse_int_after_prefix(arg, "--max-index=", options.max_index)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.max_index >= 10;
}

i64 pentagonal(const i64 n) {
    return n * (3 * n - 1) / 2;
}

bool is_pentagonal(const i64 x) {
    if (x <= 0) {
        return false;
    }
    const long double disc = 1.0L + 24.0L * static_cast<long double>(x);
    const long double root = std::sqrt(disc);
    const long double n = (1.0L + root) / 6.0L;
    const i64 k = static_cast<i64>(std::llround(n));
    return pentagonal(k) == x;
}

i64 solve(const int max_index) {
    std::vector<i64> p(static_cast<std::size_t>(max_index + 1));
    for (int i = 1; i <= max_index; ++i) {
        p[static_cast<std::size_t>(i)] = pentagonal(i);
    }

    i64 best = std::numeric_limits<i64>::max();

    for (int j = 2; j <= max_index; ++j) {
        for (int k = j - 1; k >= 1; --k) {
            const i64 diff = p[static_cast<std::size_t>(j)] - p[static_cast<std::size_t>(k)];
            if (diff >= best) {
                break;
            }
            const i64 sum = p[static_cast<std::size_t>(j)] + p[static_cast<std::size_t>(k)];
            if (is_pentagonal(diff) && is_pentagonal(sum)) {
                best = diff;
            }
        }
    }

    return best;
}

bool run_checkpoints() {
    if (!is_pentagonal(12) || !is_pentagonal(35) || is_pentagonal(36)) {
        std::cerr << "Checkpoint failed for pentagonal predicate" << '\n';
        return false;
    }
    if (!is_pentagonal(40755)) {
        std::cerr << "Checkpoint failed for known pentagonal value" << '\n';
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

    std::cout << solve(options.max_index) << '\n';
    return 0;
}
