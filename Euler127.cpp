#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 120000;
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
    return options.limit >= 3;
}

i64 solve(const int limit) {
    std::vector<int> rad(static_cast<std::size_t>(limit + 1), 1);
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    is_prime[0] = false;
    is_prime[1] = false;

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

    std::vector<std::pair<int, int>> by_rad;
    by_rad.reserve(static_cast<std::size_t>(limit));
    for (int n = 1; n < limit; ++n) {
        by_rad.push_back({rad[static_cast<std::size_t>(n)], n});
    }
    std::sort(by_rad.begin(), by_rad.end());

    i64 sum_c = 0;

    for (int c = 3; c < limit; ++c) {
        const int rc = rad[static_cast<std::size_t>(c)];
        const int threshold = c / rc;

        for (const auto& [ra, a] : by_rad) {
            if (ra >= threshold) {
                break;
            }
            if (a >= c / 2) {
                continue;
            }

            const int b = c - a;
            const i64 rad_product = static_cast<i64>(ra) *
                                    static_cast<i64>(rad[static_cast<std::size_t>(b)]) *
                                    static_cast<i64>(rc);
            if (rad_product >= c) {
                continue;
            }
            if (std::gcd(a, b) != 1) {
                continue;
            }
            sum_c += c;
        }
    }

    return sum_c;
}

bool run_checkpoints() {
    if (solve(1000) != 12523) {
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
