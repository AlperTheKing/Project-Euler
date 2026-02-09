#include <array>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int kTargetSum = 23;
constexpr int kMod23 = 23;
constexpr u64 kDefaultMod = 1000000000ULL;
constexpr u64 kDefaultN = 3138428376721ULL;  // 11^12

struct Options {
    u64 n = kDefaultN;
    u64 mod = kDefaultMod;
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.mod > 0;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

using Coeff = std::array<u64, kTargetSum + 1>;

Coeff poly_mul_mod(const Coeff& a, const Coeff& b, const u64 mod) {
    Coeff c{};
    c.fill(0);

    for (int i = 0; i <= kTargetSum; ++i) {
        if (a[static_cast<std::size_t>(i)] == 0) {
            continue;
        }
        for (int j = 0; i + j <= kTargetSum; ++j) {
            if (b[static_cast<std::size_t>(j)] == 0) {
                continue;
            }
            const u128 add = static_cast<u128>(a[static_cast<std::size_t>(i)]) *
                             static_cast<u128>(b[static_cast<std::size_t>(j)]);
            c[static_cast<std::size_t>(i + j)] =
                (c[static_cast<std::size_t>(i + j)] + static_cast<u64>(add % mod)) % mod;
        }
    }
    return c;
}

Coeff poly_pow_digit_sum(const u64 exponent, const u64 mod) {
    Coeff base{};
    base.fill(0);
    for (int d = 0; d <= 9; ++d) {
        base[static_cast<std::size_t>(d)] = 1;
    }

    Coeff result{};
    result.fill(0);
    result[0] = 1;

    u64 e = exponent;
    while (e > 0) {
        if ((e & 1ULL) != 0ULL) {
            result = poly_mul_mod(result, base, mod);
        }
        e >>= 1U;
        if (e > 0) {
            base = poly_mul_mod(base, base, mod);
        }
    }

    return result;
}

u64 count_via_classes_mod(const u64 n, const u64 mod) {
    const u64 q = n / 22ULL;
    const int r = static_cast<int>(n % 22ULL);

    const Coeff c_q = poly_pow_digit_sum(q, mod);
    const Coeff c_q1 = poly_pow_digit_sum(q + 1ULL, mod);

    std::array<int, 22> weights{};
    int cur = 1;
    for (int i = 0; i < 22; ++i) {
        weights[static_cast<std::size_t>(i)] = cur;
        cur = (cur * 10) % kMod23;
    }

    std::array<std::array<u64, kMod23>, kTargetSum + 1> dp{};
    std::array<std::array<u64, kMod23>, kTargetSum + 1> next{};
    for (auto& row : dp) {
        row.fill(0);
    }
    dp[0][0] = 1;

    for (int cls = 0; cls < 22; ++cls) {
        const Coeff& coeff = (cls < r) ? c_q1 : c_q;
        for (auto& row : next) {
            row.fill(0);
        }

        for (int s = 0; s <= kTargetSum; ++s) {
            for (int rem = 0; rem < kMod23; ++rem) {
                const u64 cur_count = dp[static_cast<std::size_t>(s)][static_cast<std::size_t>(rem)];
                if (cur_count == 0) {
                    continue;
                }
                for (int t = 0; s + t <= kTargetSum; ++t) {
                    const u64 ways = coeff[static_cast<std::size_t>(t)];
                    if (ways == 0) {
                        continue;
                    }
                    const int rem2 = (rem + weights[static_cast<std::size_t>(cls)] * t) % kMod23;
                    const u128 add = static_cast<u128>(cur_count) * static_cast<u128>(ways);
                    next[static_cast<std::size_t>(s + t)][static_cast<std::size_t>(rem2)] =
                        (next[static_cast<std::size_t>(s + t)][static_cast<std::size_t>(rem2)] +
                         static_cast<u64>(add % mod)) % mod;
                }
            }
        }

        dp = next;
    }

    return dp[kTargetSum][0];
}

u128 direct_count_exact(const int n) {
    std::array<std::array<u128, kMod23>, kTargetSum + 1> dp{};
    std::array<std::array<u128, kMod23>, kTargetSum + 1> next{};
    for (auto& row : dp) {
        row.fill(0);
    }
    dp[0][0] = 1;

    int pow10 = 1;
    for (int pos = 0; pos < n; ++pos) {
        for (auto& row : next) {
            row.fill(0);
        }

        for (int s = 0; s <= kTargetSum; ++s) {
            for (int rem = 0; rem < kMod23; ++rem) {
                const u128 cur_count = dp[static_cast<std::size_t>(s)][static_cast<std::size_t>(rem)];
                if (cur_count == 0) {
                    continue;
                }
                for (int d = 0; s + d <= kTargetSum && d <= 9; ++d) {
                    const int rem2 = (rem + d * pow10) % kMod23;
                    next[static_cast<std::size_t>(s + d)][static_cast<std::size_t>(rem2)] += cur_count;
                }
            }
        }

        dp = next;
        pow10 = (pow10 * 10) % kMod23;
    }

    return dp[kTargetSum][0];
}

u64 direct_count_mod(const int n, const u64 mod) {
    std::array<std::array<u64, kMod23>, kTargetSum + 1> dp{};
    std::array<std::array<u64, kMod23>, kTargetSum + 1> next{};
    for (auto& row : dp) {
        row.fill(0);
    }
    dp[0][0] = 1;

    int pow10 = 1;
    for (int pos = 0; pos < n; ++pos) {
        for (auto& row : next) {
            row.fill(0);
        }

        for (int s = 0; s <= kTargetSum; ++s) {
            for (int rem = 0; rem < kMod23; ++rem) {
                const u64 cur_count = dp[static_cast<std::size_t>(s)][static_cast<std::size_t>(rem)];
                if (cur_count == 0) {
                    continue;
                }
                for (int d = 0; s + d <= kTargetSum && d <= 9; ++d) {
                    const int rem2 = (rem + d * pow10) % kMod23;
                    next[static_cast<std::size_t>(s + d)][static_cast<std::size_t>(rem2)] += cur_count;
                    if (next[static_cast<std::size_t>(s + d)][static_cast<std::size_t>(rem2)] >= mod) {
                        next[static_cast<std::size_t>(s + d)][static_cast<std::size_t>(rem2)] %= mod;
                    }
                }
            }
        }

        dp = next;
        pow10 = (pow10 * 10) % kMod23;
    }

    return dp[kTargetSum][0] % mod;
}

bool run_checkpoints() {
    const u128 s9 = direct_count_exact(9);
    if (s9 != static_cast<u128>(263626ULL)) {
        std::cerr << "Checkpoint failed for S(9): got " << to_string_u128(s9) << '\n';
        return false;
    }

    const u128 s42 = direct_count_exact(42);
    if (s42 != static_cast<u128>(6377168878570056ULL)) {
        std::cerr << "Checkpoint failed for S(42): got " << to_string_u128(s42) << '\n';
        return false;
    }

    const u64 direct200 = direct_count_mod(200, kDefaultMod);
    const u64 class200 = count_via_classes_mod(200, kDefaultMod);
    if (direct200 != class200) {
        std::cerr << "Class-DP mismatch at n=200: direct=" << direct200
                  << ", class=" << class200 << '\n';
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

    const u64 answer = count_via_classes_mod(options.n, options.mod);
    std::cout << answer << '\n';
    return 0;
}
