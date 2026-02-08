#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 limit = 1000000000000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

int gcd_int(int a, int b) {
    while (b != 0) {
        const int t = a % b;
        a = b;
        b = t;
    }
    return a;
}

bool is_square(const u64 value) {
    const u64 root = static_cast<u64>(std::sqrt(static_cast<long double>(value)));
    return root * root == value || (root + 1) * (root + 1) == value;
}

u64 solve(const u64 limit) {
    std::set<u64> progressive_squares;

    for (u64 a = 1; a * a <= limit; ++a) {
        for (u64 b = 1; b < a && a * a * b * b <= limit; ++b) {
            if (gcd_int(static_cast<int>(a), static_cast<int>(b)) != 1) {
                continue;
            }

            for (u64 k = 1; k * k * a * a * b * b <= limit; ++k) {
                const u64 candidate1 = k * a * a * k * b * b + k * a * b;
                if (candidate1 <= limit && is_square(candidate1)) {
                    progressive_squares.insert(candidate1);
                }

                const u64 candidate2 = k * a * a * k * a * b + k * b * b;
                if (candidate2 <= limit && is_square(candidate2)) {
                    progressive_squares.insert(candidate2);
                }
            }
        }
    }

    u64 sum = 0;
    for (u64 value : progressive_squares) {
        sum += value;
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(100000ULL) != 124657ULL) {
        std::cerr << "Checkpoint failed for limit 100000" << '\n';
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
