#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct Options {
    int a = 25;
    int b = 75;
    int c = 1984;
    u64 mod = 100000000ULL;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(ch - '0');
        if (parsed > static_cast<i64>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    value = static_cast<int>(parsed);
    return true;
}

bool parse_u64_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--a=", options.a)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--b=", options.b)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--c=", options.c)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

std::string to_string_i128(i128 value) {
    if (value == 0) {
        return "0";
    }

    bool negative = false;
    u128 magnitude = 0;
    if (value < 0) {
        negative = true;
        magnitude = static_cast<u128>(-value);
    } else {
        magnitude = static_cast<u128>(value);
    }

    std::string digits;
    while (magnitude > 0) {
        const int digit = static_cast<int>(magnitude % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        magnitude /= 10U;
    }

    if (negative) {
        digits.push_back('-');
    }
    std::reverse(digits.begin(), digits.end());
    return digits;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) %
                            static_cast<u128>(mod));
}

u64 pow_mod(u64 base, int exp, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }

    u64 result = 1ULL % mod;
    base %= mod;

    while (exp > 0) {
        if ((exp & 1) != 0) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }

    return result;
}

i128 choose_exact(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    i128 result = 1;
    for (int i = 1; i <= k; ++i) {
        result = (result * static_cast<i128>(n - k + i)) / static_cast<i128>(i);
    }
    return result;
}

i128 pow_i128(i128 base, int exp) {
    i128 result = 1;
    while (exp > 0) {
        if ((exp & 1) != 0) {
            result *= base;
        }
        base *= base;
        exp >>= 1;
    }
    return result;
}

u64 count_unit_extensions_bruteforce(const int c, const bool is_type_a) {
    if (c < 2) {
        return 0ULL;
    }

    constexpr int left_top = 0;
    constexpr int left_bottom = 1;

    u64 count = 0ULL;

    for (int right_top = 0; right_top < c; ++right_top) {
        if (right_top == left_top) {
            continue;
        }

        for (int right_bottom = 0; right_bottom < c; ++right_bottom) {
            if (right_bottom == right_top) {
                continue;
            }
            if (is_type_a && right_bottom == left_bottom) {
                continue;
            }

            for (int upper = 0; upper < c; ++upper) {
                if (upper == left_top || upper == right_top) {
                    continue;
                }

                for (int middle = 0; middle < c; ++middle) {
                    if (middle == upper) {
                        continue;
                    }

                    for (int lower = 0; lower < c; ++lower) {
                        if (lower == middle || lower == left_bottom ||
                            lower == right_bottom) {
                            continue;
                        }
                        ++count;
                    }
                }
            }
        }
    }

    return count;
}

std::array<i128, 6> build_forward_differences(const bool is_type_a) {
    std::array<i128, 6> row{};
    for (int i = 0; i < 6; ++i) {
        row[static_cast<std::size_t>(i)] = static_cast<i128>(
            count_unit_extensions_bruteforce(i + 2, is_type_a));
    }

    std::array<i128, 6> differences{};
    differences[0] = row[0];
    for (int level = 1; level < 6; ++level) {
        for (int i = 0; i < 6 - level; ++i) {
            row[static_cast<std::size_t>(i)] =
                row[static_cast<std::size_t>(i + 1)] -
                row[static_cast<std::size_t>(i)];
        }
        differences[static_cast<std::size_t>(level)] = row[0];
    }

    return differences;
}

i128 evaluate_from_forward_differences(const std::array<i128, 6>& differences,
                                       const int c) {
    if (c < 2) {
        return 0;
    }

    const i64 n = static_cast<i64>(c) - 2LL;
    i128 choose_n_k = 1;
    i128 value = 0;

    for (int k = 0; k < 6; ++k) {
        if (k > 0) {
            choose_n_k = (choose_n_k * static_cast<i128>(n - (k - 1))) /
                         static_cast<i128>(k);
        }
        value += differences[static_cast<std::size_t>(k)] * choose_n_k;
    }

    return value;
}

i128 unit_extensions(const int c, const bool is_type_a) {
    static const std::array<i128, 6> a_differences =
        build_forward_differences(true);
    static const std::array<i128, 6> b_differences =
        build_forward_differences(false);

    return evaluate_from_forward_differences(is_type_a ? a_differences
                                                       : b_differences,
                                             c);
}

u64 normalize_mod(const i128 value, const u64 mod) {
    i128 reduced = value % static_cast<i128>(mod);
    if (reduced < 0) {
        reduced += static_cast<i128>(mod);
    }
    return static_cast<u64>(reduced);
}

i128 count_configurations_exact_small(const int a, const int b, const int c) {
    if (c < 2) {
        return 0;
    }

    const i128 a_extensions = unit_extensions(c, true);
    const i128 b_extensions = unit_extensions(c, false);

    i128 result = choose_exact(a + b, a);
    result *= static_cast<i128>(c) * static_cast<i128>(c - 1);
    result *= pow_i128(a_extensions, a);
    result *= pow_i128(b_extensions, b);
    return result;
}

u64 count_configurations_mod(const int a,
                             const int b,
                             const int c,
                             const u64 mod) {
    if (mod == 0ULL) {
        throw std::runtime_error("Modulus must be positive");
    }
    if (c < 2) {
        return 0ULL;
    }

    const i128 a_extensions = unit_extensions(c, true);
    const i128 b_extensions = unit_extensions(c, false);

    u64 result = static_cast<u64>(choose_exact(a + b, a) % static_cast<i128>(mod));

    const u64 c_mod = static_cast<u64>(static_cast<i128>(c) % static_cast<i128>(mod));
    const u64 c_minus_one_mod =
        static_cast<u64>(static_cast<i128>(c - 1) % static_cast<i128>(mod));
    result = mul_mod(result, mul_mod(c_mod, c_minus_one_mod, mod), mod);
    result = mul_mod(result, pow_mod(normalize_mod(a_extensions, mod), a, mod), mod);
    result = mul_mod(result, pow_mod(normalize_mod(b_extensions, mod), b, mod), mod);

    return result;
}

bool check_equal(const i128 actual,
                 const i128 expected,
                 const std::string& label) {
    if (actual == expected) {
        return true;
    }
    std::cerr << "Checkpoint failed for " << label << ": expected "
              << to_string_i128(expected) << ", got " << to_string_i128(actual)
              << '\n';
    return false;
}

bool run_checkpoints() {
    for (int c = 2; c <= 10; ++c) {
        const i128 a_interp = unit_extensions(c, true);
        const i128 b_interp = unit_extensions(c, false);
        const i128 a_brute =
            static_cast<i128>(count_unit_extensions_bruteforce(c, true));
        const i128 b_brute =
            static_cast<i128>(count_unit_extensions_bruteforce(c, false));

        if (!check_equal(a_interp, a_brute,
                         "unit A interpolation at c=" + std::to_string(c))) {
            return false;
        }
        if (!check_equal(b_interp, b_brute,
                         "unit B interpolation at c=" + std::to_string(c))) {
            return false;
        }
    }

    if (!check_equal(count_configurations_exact_small(1, 0, 3), 24,
                     "N(1,0,3)")) {
        return false;
    }
    if (!check_equal(count_configurations_exact_small(0, 2, 4), 92928,
                     "N(0,2,4)")) {
        return false;
    }
    if (!check_equal(count_configurations_exact_small(2, 2, 3), 20736,
                     "N(2,2,3)")) {
        return false;
    }

    return true;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }

        if (options.a < 0 || options.b < 0 || options.c < 0) {
            std::cerr << "Parameters --a, --b, and --c must be non-negative.\n";
            return 1;
        }

        if (options.run_checkpoints && !run_checkpoints()) {
            return 1;
        }

        const u64 answer =
            count_configurations_mod(options.a, options.b, options.c, options.mod);

        if (options.mod == 100000000ULL) {
            std::cout << std::setw(8) << std::setfill('0') << answer << '\n';
        } else {
            std::cout << answer << '\n';
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
