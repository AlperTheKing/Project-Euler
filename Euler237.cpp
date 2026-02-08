#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    u64 n = 1000000000000ULL;
    int modulo = 100000000;
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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

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
        if (parse_u64_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--modulo=", options.modulo)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL && options.modulo >= 2;
}

i64 mod_norm(const i64 x, const int mod) {
    i64 v = x % mod;
    if (v < 0) {
        v += mod;
    }
    return v;
}

using Mat4 = std::array<std::array<i64, 4>, 4>;
using Vec4 = std::array<i64, 4>;

Mat4 mat_mul(const Mat4& a, const Mat4& b, const int mod) {
    Mat4 c{};
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            if (a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] == 0) {
                continue;
            }
            for (int j = 0; j < 4; ++j) {
                c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                    mod_norm(c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                                 a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] *
                                     b[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)],
                             mod);
            }
        }
    }
    return c;
}

Vec4 mat_vec_mul(const Mat4& a, const Vec4& x, const int mod) {
    Vec4 y{};
    for (int i = 0; i < 4; ++i) {
        i64 sum = 0;
        for (int j = 0; j < 4; ++j) {
            sum += a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
                   x[static_cast<std::size_t>(j)];
        }
        y[static_cast<std::size_t>(i)] = mod_norm(sum, mod);
    }
    return y;
}

i64 solve(const u64 n, const int mod) {
    // T(1)=1, T(2)=1, T(3)=4, T(4)=8
    if (n == 1ULL || n == 2ULL) {
        return 1 % mod;
    }
    if (n == 3ULL) {
        return 4 % mod;
    }
    if (n == 4ULL) {
        return 8 % mod;
    }

    Mat4 step{};
    step[0] = {2, 2, mod_norm(-2, mod), 1};
    step[1] = {1, 0, 0, 0};
    step[2] = {0, 1, 0, 0};
    step[3] = {0, 0, 1, 0};

    Mat4 acc{};
    for (int i = 0; i < 4; ++i) {
        acc[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1;
    }

    u64 power = n - 4ULL;
    while (power > 0ULL) {
        if ((power & 1ULL) != 0ULL) {
            acc = mat_mul(acc, step, mod);
        }
        step = mat_mul(step, step, mod);
        power >>= 1U;
    }

    const Vec4 base{8, 4, 1, 1};  // [T4, T3, T2, T1]
    const Vec4 out = mat_vec_mul(acc, base, mod);
    return out[0];
}

u64 brute_tours(int n) {
    const int rows = 4;
    const int cols = n;
    const int total = rows * cols;
    if (total > 20) {
        return 0;
    }

    const auto vid = [cols](const int r, const int c) { return r * cols + c; };
    std::vector<std::vector<int>> g(static_cast<std::size_t>(total));
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const int v = vid(r, c);
            if (r > 0) {
                g[static_cast<std::size_t>(v)].push_back(vid(r - 1, c));
            }
            if (r + 1 < rows) {
                g[static_cast<std::size_t>(v)].push_back(vid(r + 1, c));
            }
            if (c > 0) {
                g[static_cast<std::size_t>(v)].push_back(vid(r, c - 1));
            }
            if (c + 1 < cols) {
                g[static_cast<std::size_t>(v)].push_back(vid(r, c + 1));
            }
        }
    }

    const int start = vid(0, 0);
    const int target = vid(3, 0);
    const int max_mask = 1 << total;
    std::vector<std::vector<u64>> dp(static_cast<std::size_t>(max_mask),
                                     std::vector<u64>(static_cast<std::size_t>(total), 0ULL));
    dp[1 << start][static_cast<std::size_t>(start)] = 1ULL;

    for (int mask = 0; mask < max_mask; ++mask) {
        for (int v = 0; v < total; ++v) {
            const u64 ways = dp[static_cast<std::size_t>(mask)][static_cast<std::size_t>(v)];
            if (ways == 0) {
                continue;
            }
            for (int to : g[static_cast<std::size_t>(v)]) {
                if ((mask & (1 << to)) != 0) {
                    continue;
                }
                const int nmask = mask | (1 << to);
                dp[static_cast<std::size_t>(nmask)][static_cast<std::size_t>(to)] += ways;
            }
        }
    }
    return dp[static_cast<std::size_t>(max_mask - 1)][static_cast<std::size_t>(target)];
}

bool run_checkpoints() {
    for (int n = 1; n <= 4; ++n) {
        if (solve(static_cast<u64>(n), 1000000000) != static_cast<i64>(brute_tours(n))) {
            std::cerr << "Checkpoint failed for brute cross-check at n=" << n << '\n';
            return false;
        }
    }
    if (solve(10, 1000000000) != 2329) {
        std::cerr << "Checkpoint failed for T(10)=2329 sample" << '\n';
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
    std::cout << solve(options.n, options.modulo) << '\n';
    return 0;
}
