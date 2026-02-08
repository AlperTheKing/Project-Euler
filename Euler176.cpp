#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int target_count = 47547;
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
        if (parse_int_after_prefix(arg, "--target-count=", options.target_count)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target_count >= 1;
}

u64 pow_u64(u64 base, int exp) {
    u64 result = 1;
    while (exp-- > 0) {
        result *= base;
    }
    return result;
}

u64 divisor_count_square(u64 n) {
    u64 result = 1;
    int exponent = 0;
    while ((n & 1ULL) == 0ULL) {
        n >>= 1U;
        ++exponent;
    }
    if (exponent > 0) {
        result *= static_cast<u64>(2 * exponent + 1);
    }
    for (u64 p = 3; p * p <= n; p += 2) {
        if (n % p != 0) {
            continue;
        }
        exponent = 0;
        while (n % p == 0) {
            n /= p;
            ++exponent;
        }
        result *= static_cast<u64>(2 * exponent + 1);
    }
    if (n > 1) {
        result *= 3ULL;
    }
    return result;
}

u64 count_triangles_for_cathetus(const u64 n) {
    const u64 m = (n % 2ULL == 0ULL) ? (n / 2ULL) : n;
    return (divisor_count_square(m) - 1ULL) / 2ULL;
}

void collect_factorizations(int remaining, int min_factor,
                            std::vector<int>& current,
                            std::vector<std::vector<int>>& out) {
    if (remaining == 1) {
        out.push_back(current);
        return;
    }
    for (int factor = min_factor; factor <= remaining; factor += 2) {
        if (remaining % factor != 0) {
            continue;
        }
        current.push_back(factor);
        collect_factorizations(remaining / factor, factor, current, out);
        current.pop_back();
    }
}

u128 build_number(int two_exponent, std::vector<int> exponents) {
    static const std::vector<u64> odd_primes{
        3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL,
        37ULL, 41ULL, 43ULL, 47ULL, 53ULL, 59ULL, 61ULL, 67ULL, 71ULL, 73ULL};
    std::sort(exponents.begin(), exponents.end(), std::greater<int>());

    u128 value = 1;
    if (two_exponent > 0) {
        value <<= two_exponent;
    }
    for (std::size_t i = 0; i < exponents.size(); ++i) {
        value *= static_cast<u128>(pow_u64(odd_primes[i], exponents[i]));
    }
    return value;
}

u128 solve_u128(const int target_count) {
    const int product_target = 2 * target_count + 1;
    u128 best = std::numeric_limits<u128>::max();

    auto try_odd_factorization = [&](int two_factor_odd_part) {
        if (product_target % two_factor_odd_part != 0) {
            return;
        }
        const int remaining = product_target / two_factor_odd_part;
        const int two_exponent = (two_factor_odd_part + 1) / 2;

        std::vector<std::vector<int>> factorizations;
        std::vector<int> current;
        collect_factorizations(remaining, 3, current, factorizations);
        if (remaining == 1) {
            factorizations.push_back({});
        }

        for (const auto& odd_factors : factorizations) {
            std::vector<int> odd_exponents;
            odd_exponents.reserve(odd_factors.size());
            for (int f : odd_factors) {
                odd_exponents.push_back((f - 1) / 2);
            }
            const u128 candidate = build_number(two_exponent, odd_exponents);
            if (candidate < best) {
                best = candidate;
            }
        }
    };

    // Odd catheti: no factor from powers of 2 in the divisor-count product.
    {
        std::vector<std::vector<int>> factorizations;
        std::vector<int> current;
        collect_factorizations(product_target, 3, current, factorizations);
        for (const auto& odd_factors : factorizations) {
            std::vector<int> odd_exponents;
            odd_exponents.reserve(odd_factors.size());
            for (int f : odd_factors) {
                odd_exponents.push_back((f - 1) / 2);
            }
            const u128 candidate = build_number(0, odd_exponents);
            if (candidate < best) {
                best = candidate;
            }
        }
    }

    // Even catheti: (2k-1) contributes to the divisor-count product.
    for (int factor = 1; factor <= product_target; factor += 2) {
        if (product_target % factor != 0) {
            continue;
        }
        try_odd_factorization(factor);
    }

    return best;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

bool run_checkpoints() {
    if (count_triangles_for_cathetus(12) != 4ULL) {
        std::cerr << "Checkpoint failed for cathetus 12" << '\n';
        return false;
    }
    if (to_string_u128(solve_u128(4)) != "12") {
        std::cerr << "Checkpoint failed for target-count 4" << '\n';
        return false;
    }
    if (to_string_u128(solve_u128(1)) != "3") {
        std::cerr << "Checkpoint failed for target-count 1" << '\n';
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
    std::cout << to_string_u128(solve_u128(options.target_count)) << '\n';
    return 0;
}
