#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <tuple>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1'000'000'123ULL;

struct Options {
    int max_power = 39;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
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
        if (parse_int_after_prefix(arg, "--max-power=", options.max_power)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_power >= 1 && options.max_power <= 62;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    base %= MOD;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * base) % MOD);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % MOD);
        exp >>= 1ULL;
    }
    return result;
}

u64 mod_inverse(const u64 x) {
    return mod_pow(x % MOD, MOD - 2ULL);
}

u64 geom_sum_from_one(const i64 lambda, const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    if (lambda == 1) {
        return n % MOD;
    }
    if (lambda == 0) {
        return 0ULL;
    }
    const u64 lm = static_cast<u64>((lambda % static_cast<i64>(MOD) + static_cast<i64>(MOD)) % static_cast<i64>(MOD));
    const u64 num = static_cast<u64>((static_cast<u128>(lm) * ((mod_pow(lm, n) + MOD - 1ULL) % MOD)) % MOD);
    const u64 den = (lm + MOD - 1ULL) % MOD;
    return static_cast<u64>((static_cast<u128>(num) * mod_inverse(den)) % MOD);
}

u64 geom_sum_from_zero(const i64 lambda, const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    if (lambda == 1) {
        return n % MOD;
    }
    if (lambda == 0) {
        return 1ULL;
    }
    const u64 lm = static_cast<u64>((lambda % static_cast<i64>(MOD) + static_cast<i64>(MOD)) % static_cast<i64>(MOD));
    const u64 num = (mod_pow(lm, n) + MOD - 1ULL) % MOD;
    const u64 den = (lm + MOD - 1ULL) % MOD;
    return static_cast<u64>((static_cast<u128>(num) * mod_inverse(den)) % MOD);
}

std::map<i64, i64> expand_coeffs_F() {
    std::map<std::pair<int, int>, i64> poly1;
    poly1[{0, 0}] = 1;
    const std::array<std::tuple<int, int, i64>, 3> base1 = {
        std::make_tuple(0, 0, 2), std::make_tuple(1, 0, 1), std::make_tuple(0, 1, -1)};

    for (int rep = 0; rep < 5; ++rep) {
        std::map<std::pair<int, int>, i64> next;
        for (const auto& row : poly1) {
            for (const auto& term : base1) {
                const int da = std::get<0>(term);
                const int db = std::get<1>(term);
                const i64 dc = std::get<2>(term);
                next[{row.first.first + da, row.first.second + db}] += row.second * dc;
            }
        }
        poly1.swap(next);
    }

    std::map<std::pair<int, int>, i64> poly2;
    poly2[{0, 0}] = 1;
    const std::array<std::tuple<int, int, i64>, 2> base2 = {std::make_tuple(1, 0, 1), std::make_tuple(0, 1, 1)};
    for (int rep = 0; rep < 5; ++rep) {
        std::map<std::pair<int, int>, i64> next;
        for (const auto& row : poly2) {
            for (const auto& term : base2) {
                next[{row.first.first + std::get<0>(term), row.first.second + std::get<1>(term)}] +=
                    row.second * std::get<2>(term);
            }
        }
        poly2.swap(next);
    }

    std::map<std::pair<int, int>, i64> both;
    for (const auto& a : poly1) {
        for (const auto& b : poly2) {
            both[{a.first.first + b.first.first, a.first.second + b.first.second}] += a.second * b.second;
        }
    }

    std::map<i64, i64> by_lambda;
    for (const auto& row : both) {
        by_lambda[static_cast<i64>(row.first.first - row.first.second)] += row.second;
    }
    return by_lambda;
}

std::map<i64, i64> expand_coeffs_G() {
    std::map<std::pair<int, int>, i64> poly1;
    poly1[{0, 0}] = 1;
    const std::array<std::tuple<int, int, i64>, 3> base1 = {
        std::make_tuple(0, 0, 2), std::make_tuple(1, 0, 1), std::make_tuple(0, 1, -1)};
    for (int rep = 0; rep < 5; ++rep) {
        std::map<std::pair<int, int>, i64> next;
        for (const auto& row : poly1) {
            for (const auto& term : base1) {
                next[{row.first.first + std::get<0>(term), row.first.second + std::get<1>(term)}] +=
                    row.second * std::get<2>(term);
            }
        }
        poly1.swap(next);
    }

    std::map<std::pair<int, int>, i64> poly2;
    poly2[{0, 0}] = 1;
    const std::array<std::tuple<int, int, i64>, 2> plusminus = {std::make_tuple(1, 0, 1), std::make_tuple(0, 1, 1)};
    for (int rep = 0; rep < 4; ++rep) {
        std::map<std::pair<int, int>, i64> next;
        for (const auto& row : poly2) {
            for (const auto& term : plusminus) {
                next[{row.first.first + std::get<0>(term), row.first.second + std::get<1>(term)}] +=
                    row.second * std::get<2>(term);
            }
        }
        poly2.swap(next);
    }

    std::map<std::pair<int, int>, i64> diff;
    diff[{1, 0}] = 1;
    diff[{0, 1}] = -1;

    std::map<std::pair<int, int>, i64> tmp;
    for (const auto& a : poly1) {
        for (const auto& b : poly2) {
            tmp[{a.first.first + b.first.first, a.first.second + b.first.second}] += a.second * b.second;
        }
    }
    std::map<std::pair<int, int>, i64> both;
    for (const auto& a : tmp) {
        for (const auto& b : diff) {
            both[{a.first.first + b.first.first, a.first.second + b.first.second}] += a.second * b.second;
        }
    }

    std::map<i64, i64> by_lambda;
    for (const auto& row : both) {
        by_lambda[static_cast<i64>(row.first.first - row.first.second)] += row.second;
    }
    return by_lambda;
}

u64 Q(const u64 n, const std::map<i64, i64>& cf, const std::map<i64, i64>& cg, const u64 inv_two_pow_10) {
    u64 sum_A = 0ULL;
    for (const auto& row : cf) {
        const i64 coeff = row.second;
        const u64 geom = geom_sum_from_one(row.first, n);
        const u64 coeff_mod = static_cast<u64>((coeff % static_cast<i64>(MOD) + static_cast<i64>(MOD)) %
                                               static_cast<i64>(MOD));
        sum_A = (sum_A + static_cast<u64>((static_cast<u128>(coeff_mod) * geom) % MOD)) % MOD;
    }

    u64 sum_B = 0ULL;
    for (const auto& row : cg) {
        const i64 coeff = row.second;
        const u64 geom = geom_sum_from_zero(row.first, n);
        const u64 coeff_mod = static_cast<u64>((coeff % static_cast<i64>(MOD) + static_cast<i64>(MOD)) %
                                               static_cast<i64>(MOD));
        sum_B = (sum_B + static_cast<u64>((static_cast<u128>(coeff_mod) * geom) % MOD)) % MOD;
    }

    return static_cast<u64>((static_cast<u128>((sum_A + MOD - sum_B) % MOD) * inv_two_pow_10) % MOD);
}

u64 solve(const int max_power) {
    const std::map<i64, i64> cf = expand_coeffs_F();
    const std::map<i64, i64> cg = expand_coeffs_G();
    const u64 inv_two_pow_10 = mod_inverse(1024ULL);

    u64 ans = 0ULL;
    for (int u = 1; u <= max_power; ++u) {
        const u64 n = 1ULL << static_cast<unsigned>(u);
        ans += Q(n, cf, cg, inv_two_pow_10);
        if (ans >= MOD) {
            ans -= MOD;
        }
    }
    return ans;
}

u64 brute_Q_small(const int n) {
    u64 count = 0ULL;
    const auto is_simber = [](u64 x) {
        std::array<int, 10> cnt{};
        while (x > 0ULL) {
            ++cnt[static_cast<std::size_t>(x % 10ULL)];
            x /= 10ULL;
        }
        for (int d = 0; d <= 9; ++d) {
            if ((d % 2) == 0) {
                if ((cnt[static_cast<std::size_t>(d)] % 2) != 0) {
                    return false;
                }
            } else {
                if ((cnt[static_cast<std::size_t>(d)] % 2) == 0 && cnt[static_cast<std::size_t>(d)] != 0) {
                    return false;
                }
            }
        }
        return true;
    };

    u64 upper = 1ULL;
    for (int i = 0; i < n; ++i) {
        upper *= 10ULL;
    }
    for (u64 x = 1ULL; x < upper; ++x) {
        if (is_simber(x)) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    const std::map<i64, i64> cf = expand_coeffs_F();
    const std::map<i64, i64> cg = expand_coeffs_G();
    const u64 inv_two_pow_10 = mod_inverse(1024ULL);

    if (Q(7ULL, cf, cg, inv_two_pow_10) != 287'975ULL) {
        std::cerr << "Checkpoint failed: Q(7)=287975" << '\n';
        return false;
    }
    if (Q(100ULL, cf, cg, inv_two_pow_10) != 123'864'868ULL) {
        std::cerr << "Checkpoint failed: Q(100) mod 1,000,000,123 = 123864868" << '\n';
        return false;
    }
    if (Q(4ULL, cf, cg, inv_two_pow_10) != brute_Q_small(4) % MOD) {
        std::cerr << "Checkpoint failed: brute-force cross-check for Q(4)" << '\n';
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
    std::cout << solve(options.max_power) << '\n';
    return 0;
}
