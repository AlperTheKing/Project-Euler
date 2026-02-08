#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 10000;
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

    i64 parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
        if (parsed > static_cast<i64>(std::numeric_limits<int>::max())) {
            return false;
        }
    }
    value = static_cast<int>(parsed);
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

i64 proper_divisor_sum_via_factorization(i64 n) {
    if (n <= 1) {
        return 0;
    }

    i64 remaining = n;
    i64 sigma = 1;

    int exponent = 0;
    while ((remaining & 1LL) == 0LL) {
        remaining >>= 1LL;
        ++exponent;
    }
    if (exponent > 0) {
        i64 term = 1;
        i64 pow = 1;
        for (int i = 0; i < exponent; ++i) {
            pow <<= 1LL;
            term += pow;
        }
        sigma *= term;
    }

    for (i64 p = 3; p <= remaining / p; p += 2) {
        exponent = 0;
        while (remaining % p == 0) {
            remaining /= p;
            ++exponent;
        }
        if (exponent > 0) {
            i64 term = 1;
            i64 pow = 1;
            for (int i = 0; i < exponent; ++i) {
                pow *= p;
                term += pow;
            }
            sigma *= term;
        }
    }

    if (remaining > 1) {
        sigma *= (1 + remaining);
    }

    return sigma - n;
}

std::vector<i64> proper_divisor_sums_up_to(const int limit) {
    std::vector<i64> sums(static_cast<std::size_t>(limit), 0LL);
    for (int d = 1; d <= (limit - 1) / 2; ++d) {
        for (int m = d * 2; m < limit; m += d) {
            sums[static_cast<std::size_t>(m)] += d;
        }
    }
    return sums;
}

i64 solve(const int limit) {
    const std::vector<i64> sums = proper_divisor_sums_up_to(limit);

    i64 total = 0;
    for (int a = 2; a < limit; ++a) {
        const i64 b = sums[static_cast<std::size_t>(a)];
        if (b == a || b <= 0) {
            continue;
        }

        i64 back_sum = 0;
        if (b < limit) {
            back_sum = sums[static_cast<std::size_t>(b)];
        } else {
            back_sum = proper_divisor_sum_via_factorization(b);
        }

        if (back_sum == a) {
            total += a;
        }
    }

    return total;
}

bool run_checkpoints() {
    if (proper_divisor_sum_via_factorization(220) != 284 ||
        proper_divisor_sum_via_factorization(284) != 220) {
        std::cerr << "Checkpoint failed for amicable pair (220, 284)" << '\n';
        return false;
    }
    if (solve(300) != 504) {
        std::cerr << "Checkpoint failed for limit=300" << '\n';
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
