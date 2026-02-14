#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <cmath>
#include <functional>

namespace {

using i64 = std::int64_t;

struct Options {
    i64 limit = 10000000000000000LL;
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
    return options.limit >= 1;
}

i64 brute_non_superperfect(const i64 limit) {
    i64 count = 0;

    for (i64 m = 2; m * m <= limit; ++m) {
        for (i64 n = 1; n < m; ++n) {
            if (((m - n) & 1LL) == 0LL || std::gcd(m, n) != 1LL) {
                continue;
            }

            const i64 c = m * m + n * n;
            if (c > limit) {
                break;
            }

            const i64 root = static_cast<i64>(std::sqrt(static_cast<long double>(c)));
            if (root * root != c) {
                continue;
            }

            const i64 a = m * m - n * n;
            const i64 b = 2 * m * n;
            const i64 area = (a * b) / 2;
            if (area % 84LL != 0LL) {
                ++count;
            }
        }
    }

    return count;
}

i64 solve(const i64 /*limit*/) {
    // Every primitive right triangle with square hypotenuse is super-perfect,
    // so the required count is zero.
    return 0;
}

bool run_checkpoints() {
    if (brute_non_superperfect(1000000LL) != 0LL) {
        std::cerr << "Checkpoint failed for brute bound 1e6" << '\n';
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
