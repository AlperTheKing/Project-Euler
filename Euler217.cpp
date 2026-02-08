#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int n = 47;
    u64 mod = 14348907ULL;  // 3^15
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
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.n >= 1 && options.mod >= 2;
}

u64 add_mod(const u64 a, const u64 b, const u64 mod) {
    u64 s = a + b;
    if (s >= mod) {
        s -= mod;
    }
    return s;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

struct DigitDP {
    std::vector<std::vector<u64>> count;
    std::vector<std::vector<u64>> sum;
};

DigitDP build_dp(const int max_len, const bool leading_nonzero, const u64 mod) {
    const int max_sum = 9 * max_len;
    std::vector<std::vector<u64>> count(max_len + 1, std::vector<u64>(max_sum + 1, 0));
    std::vector<std::vector<u64>> sum(max_len + 1, std::vector<u64>(max_sum + 1, 0));

    count[0][0] = 1;

    for (int len = 0; len < max_len; ++len) {
        const int cur_max_sum = 9 * len;
        for (int s = 0; s <= cur_max_sum; ++s) {
            const u64 cnt = count[len][s];
            const u64 sm = sum[len][s];
            if (cnt == 0 && sm == 0) {
                continue;
            }

            int d_start = 0;
            int d_end = 9;
            if (leading_nonzero && len == 0) {
                d_start = 1;
            }

            for (int d = d_start; d <= d_end; ++d) {
                const int ns = s + d;
                count[len + 1][ns] = add_mod(count[len + 1][ns], cnt, mod);

                u64 val = mul_mod(sm, 10ULL, mod);
                val = add_mod(val, mul_mod(static_cast<u64>(d), cnt, mod), mod);
                sum[len + 1][ns] = add_mod(sum[len + 1][ns], val, mod);
            }
        }
    }

    return {std::move(count), std::move(sum)};
}

u64 pow10_mod(const int exp, const u64 mod) {
    u64 result = 1 % mod;
    u64 base = 10 % mod;
    int e = exp;
    while (e > 0) {
        if (e & 1) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        e >>= 1;
    }
    return result;
}

u64 solve_mod(const int n, const u64 mod) {
    const int max_half = (n + 1) / 2;
    const DigitDP left_dp = build_dp(max_half, true, mod);
    const DigitDP right_dp = build_dp(max_half, false, mod);

    std::vector<u64> pow10(static_cast<std::size_t>(n + 2), 1 % mod);
    for (int i = 1; i <= n + 1; ++i) {
        pow10[static_cast<std::size_t>(i)] = mul_mod(pow10[static_cast<std::size_t>(i - 1)], 10ULL, mod);
    }

    u64 total = 0;

    for (int len = 1; len <= n; ++len) {
        const int half = len / 2;
        const bool odd = (len % 2 == 1);
        const int max_sum = 9 * half;

        const auto& cnt_left = left_dp.count[half];
        const auto& sum_left = left_dp.sum[half];
        const auto& cnt_right = right_dp.count[half];
        const auto& sum_right = right_dp.sum[half];

        if (!odd) {
            const u64 shift = pow10[static_cast<std::size_t>(half)];
            for (int s = 0; s <= max_sum; ++s) {
                const u64 cl = cnt_left[static_cast<std::size_t>(s)];
                const u64 cr = cnt_right[static_cast<std::size_t>(s)];
                if (cl == 0 || cr == 0) {
                    continue;
                }

                const u64 part_left = mul_mod(mul_mod(sum_left[static_cast<std::size_t>(s)], shift, mod), cr, mod);
                const u64 part_right = mul_mod(sum_right[static_cast<std::size_t>(s)], cl, mod);
                total = add_mod(total, add_mod(part_left, part_right, mod), mod);
            }
        } else {
            const u64 shift_left = pow10[static_cast<std::size_t>(half + 1)];
            const u64 shift_mid = pow10[static_cast<std::size_t>(half)];
            for (int s = 0; s <= max_sum; ++s) {
                const u64 cl = cnt_left[static_cast<std::size_t>(s)];
                const u64 cr = cnt_right[static_cast<std::size_t>(s)];
                if (cl == 0 || cr == 0) {
                    continue;
                }

                const u64 part_left = mul_mod(mul_mod(sum_left[static_cast<std::size_t>(s)], shift_left, mod), mul_mod(cr, 10ULL, mod), mod);
                const u64 part_mid = mul_mod(mul_mod(45ULL, shift_mid, mod), mul_mod(cl, cr, mod), mod);
                const u64 part_right = mul_mod(mul_mod(sum_right[static_cast<std::size_t>(s)], cl, mod), 10ULL, mod);

                total = add_mod(total, add_mod(add_mod(part_left, part_mid, mod), part_right, mod), mod);
            }
        }
    }

    return total;
}

u64 brute_small(const int n) {
    const u64 limit = pow10_mod(n, 1000000000000000000ULL);

    auto balanced = [](u64 x) {
        std::string s = std::to_string(x);
        const int k = static_cast<int>(s.size());
        const int m = (k + 1) / 2;
        int left = 0;
        int right = 0;
        for (int i = 0; i < m; ++i) {
            left += s[static_cast<std::size_t>(i)] - '0';
            right += s[static_cast<std::size_t>(k - 1 - i)] - '0';
        }
        return left == right;
    };

    u64 sum = 0;
    for (u64 x = 1; x < limit; ++x) {
        if (balanced(x)) {
            sum += x;
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve_mod(1, 1000000000000000000ULL) != 45ULL) {
        std::cerr << "Checkpoint failed for T(1)" << '\n';
        return false;
    }
    if (solve_mod(2, 1000000000000000000ULL) != 540ULL) {
        std::cerr << "Checkpoint failed for T(2)" << '\n';
        return false;
    }
    if (solve_mod(5, 1000000000000000000ULL) != 334795890ULL) {
        std::cerr << "Checkpoint failed for T(5)" << '\n';
        return false;
    }
    if (solve_mod(6, 1000000000000000000ULL) != brute_small(6)) {
        std::cerr << "Checkpoint failed for brute-force cross-check T(6)" << '\n';
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

    std::cout << solve_mod(options.n, options.mod) << '\n';
    return 0;
}
