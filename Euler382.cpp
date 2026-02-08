#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

constexpr u64 kDefaultN = 1'000'000'000'000'000'000ULL;
constexpr u64 kMod = 1'000'000'000ULL;
constexpr int kOrder = 8;
constexpr std::array<i64, kOrder> kRecurrence = {3, -2, 2, -5, 1, 1, 3, -2};
constexpr int kValidationN = 60;

struct Options {
    u64 n = kDefaultN;
    bool run_checkpoints = true;
};

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
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

u64 mod_norm(const i128 value, const u64 mod) {
    i128 x = value % static_cast<i128>(mod);
    if (x < 0) {
        x += static_cast<i128>(mod);
    }
    return static_cast<u64>(x);
}

u64 add_mod(const u64 a, const u64 b, const u64 mod) {
    const u64 sum = a + b;
    if (sum >= mod || sum < a) {
        return sum % mod;
    }
    return sum;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) %
                            static_cast<u128>(mod));
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }

    u64 result = 1ULL % mod;
    base %= mod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1U;
    }
    return result;
}

struct Key {
    int i = 0;
    u64 limit = 0ULL;

    bool operator==(const Key& other) const {
        return i == other.i && limit == other.limit;
    }
};

struct KeyHash {
    std::size_t operator()(const Key& key) const {
        const std::size_t h1 = std::hash<int>{}(key.i);
        const std::size_t h2 = std::hash<u64>{}(key.limit);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
    }
};

class ExactCounter {
public:
    explicit ExactCounter(const int max_n) : s_(static_cast<std::size_t>(max_n + 1), 0ULL),
                                             prefix_(static_cast<std::size_t>(max_n + 1), 0U) {
        if (max_n < 3) {
            throw std::runtime_error("ExactCounter requires max_n >= 3.");
        }

        s_[1] = 1ULL;
        s_[2] = 2ULL;
        s_[3] = 3ULL;
        for (int n = 4; n <= max_n; ++n) {
            s_[n] = s_[n - 1] + s_[n - 3];
        }
        for (int i = 1; i <= max_n; ++i) {
            prefix_[i] = prefix_[i - 1] + static_cast<u128>(s_[i]);
        }
    }

    const std::vector<u64>& sequence() const {
        return s_;
    }

    u128 count_subsets_leq(const int i, const u64 limit) {
        if (i <= 0) {
            return 1U;
        }

        if (static_cast<u128>(limit) >= prefix_[i]) {
            return static_cast<u128>(1U) << i;
        }

        const Key key{i, limit};
        const auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        u128 result = 0U;
        if (s_[i] > limit) {
            result = count_subsets_leq(i - 1, limit);
        } else {
            result = count_subsets_leq(i - 1, limit) + count_subsets_leq(i - 1, limit - s_[i]);
        }

        memo_.emplace(key, result);
        return result;
    }

private:
    std::vector<u64> s_;
    std::vector<u128> prefix_;
    std::unordered_map<Key, u128, KeyHash> memo_;
};

struct Mat9 {
    std::array<std::array<u64, 9>, 9> a{};
};

Mat9 identity9() {
    Mat9 id;
    for (int i = 0; i < 9; ++i) {
        id.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1ULL;
    }
    return id;
}

Mat9 mat_mul(const Mat9& lhs, const Mat9& rhs, const u64 mod) {
    Mat9 out;
    for (int i = 0; i < 9; ++i) {
        for (int k = 0; k < 9; ++k) {
            const u64 lv = lhs.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            if (lv == 0ULL) {
                continue;
            }
            for (int j = 0; j < 9; ++j) {
                const u64 rv = rhs.a[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
                if (rv == 0ULL) {
                    continue;
                }
                const u64 add = mul_mod(lv, rv, mod);
                u64& cell = out.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                cell = add_mod(cell, add, mod);
            }
        }
    }
    return out;
}

std::array<u64, 9> mat_vec_mul(const Mat9& m, const std::array<u64, 9>& v, const u64 mod) {
    std::array<u64, 9> out{};
    for (int i = 0; i < 9; ++i) {
        u64 acc = 0ULL;
        for (int j = 0; j < 9; ++j) {
            const u64 prod = mul_mod(m.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
                                     v[static_cast<std::size_t>(j)],
                                     mod);
            acc = add_mod(acc, prod, mod);
        }
        out[static_cast<std::size_t>(i)] = acc;
    }
    return out;
}

std::array<u64, 9> apply_mat_pow(Mat9 base,
                                 u64 exp,
                                 std::array<u64, 9> vec,
                                 const u64 mod) {
    Mat9 acc = identity9();
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            acc = mat_mul(base, acc, mod);
        }
        base = mat_mul(base, base, mod);
        exp >>= 1U;
    }
    return mat_vec_mul(acc, vec, mod);
}

u64 prefix_sum_b_mod(const u64 n, const u64 mod) {
    constexpr std::array<u64, kOrder> kBase = {1ULL, 2ULL, 4ULL, 6ULL, 11ULL, 20ULL, 36ULL,
                                                67ULL};

    if (n == 0ULL) {
        return 0ULL;
    }

    if (n <= static_cast<u64>(kOrder)) {
        u64 sum = 0ULL;
        for (u64 i = 0ULL; i < n; ++i) {
            sum = add_mod(sum, kBase[static_cast<std::size_t>(i)] % mod, mod);
        }
        return sum;
    }

    u64 s8 = 0ULL;
    for (const u64 v : kBase) {
        s8 = add_mod(s8, v % mod, mod);
    }

    std::array<u64, 9> state{
        kBase[7] % mod, kBase[6] % mod, kBase[5] % mod, kBase[4] % mod, kBase[3] % mod,
        kBase[2] % mod, kBase[1] % mod, kBase[0] % mod, s8};

    Mat9 trans;
    for (int j = 0; j < kOrder; ++j) {
        trans.a[0][static_cast<std::size_t>(j)] = mod_norm(static_cast<i128>(kRecurrence[j]), mod);
    }
    for (int r = 1; r < kOrder; ++r) {
        trans.a[static_cast<std::size_t>(r)][static_cast<std::size_t>(r - 1)] = 1ULL;
    }
    for (int j = 0; j < kOrder; ++j) {
        trans.a[8][static_cast<std::size_t>(j)] = mod_norm(static_cast<i128>(kRecurrence[j]), mod);
    }
    trans.a[8][8] = 1ULL;

    const std::array<u64, 9> advanced = apply_mat_pow(trans, n - 8ULL, state, mod);
    return advanced[8];
}

u64 solve_mod(const u64 n, const u64 mod) {
    const u64 sum_b = prefix_sum_b_mod(n, mod);
    const u64 pow2 = pow_mod(2ULL, n, mod);
    i128 ans = static_cast<i128>(pow2) - 1 - static_cast<i128>(sum_b);
    return mod_norm(ans, mod);
}

u64 to_u64(const u128 x) {
    return static_cast<u64>(x);
}

void run_checkpoints() {
    ExactCounter exact(kValidationN);
    const auto& s = exact.sequence();

    std::vector<u128> b_exact(static_cast<std::size_t>(kValidationN + 1), 0U);
    for (int n = 1; n <= kValidationN; ++n) {
        b_exact[static_cast<std::size_t>(n)] = exact.count_subsets_leq(n - 1, s[n]);
    }

    for (int n = 9; n <= kValidationN; ++n) {
        i128 rhs = 0;
        for (int j = 0; j < kOrder; ++j) {
            rhs += static_cast<i128>(kRecurrence[static_cast<std::size_t>(j)]) *
                   static_cast<i128>(b_exact[static_cast<std::size_t>(n - 1 - j)]);
        }
        const u128 lhs = b_exact[static_cast<std::size_t>(n)];
        if (rhs < 0 || static_cast<u128>(rhs) != lhs) {
            throw std::runtime_error("Derived linear recurrence failed validation.");
        }
    }

    const auto f_exact = [&](const int n) -> u128 {
        u128 sum_b = 0U;
        for (int i = 1; i <= n; ++i) {
            sum_b += b_exact[static_cast<std::size_t>(i)];
        }
        return (static_cast<u128>(1U) << n) - 1U - sum_b;
    };

    const u64 f5 = to_u64(f_exact(5));
    const u64 f10 = to_u64(f_exact(10));
    const u64 f25 = to_u64(f_exact(25));
    if (f5 != 7ULL || f10 != 501ULL || f25 != 18'635'853ULL) {
        throw std::runtime_error("Problem statement checkpoints failed.");
    }

    for (int n = 1; n <= kValidationN; ++n) {
        const u64 fast = solve_mod(static_cast<u64>(n), kMod);
        const u64 slow = static_cast<u64>(f_exact(n) % static_cast<u128>(kMod));
        if (fast != slow) {
            throw std::runtime_error("Fast modular solver mismatch on validation range.");
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    try {
        if (options.run_checkpoints) {
            run_checkpoints();
        }
    } catch (const std::exception& ex) {
        std::cerr << "Checkpoint failure: " << ex.what() << '\n';
        return 1;
    }

    const u64 answer = solve_mod(options.n, kMod);
    std::cout << "f(" << options.n << ") mod 10^9 = " << answer << '\n';
    std::cout << "Last 9 digits: " << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}
