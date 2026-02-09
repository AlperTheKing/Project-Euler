#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1307674368000ULL;  // 15!

struct Options {
    u64 n = 1000000000000000ULL;
    int x_max = 100;
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
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
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
    try {
        value = static_cast<u64>(std::stoull(tail));
    } catch (...) {
        return false;
    }
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
            parse_int_after_prefix(arg, "--x-max=", options.x_max)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.x_max >= 0;
}

using Mat3 = std::array<std::array<u64, 3>, 3>;
using Vec3 = std::array<u64, 3>;

Mat3 mat_mul(const Mat3& a, const Mat3& b) {
    Mat3 c{};
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (a[i][k] == 0ULL) {
                continue;
            }
            const u64 aik = a[i][k];
            for (int j = 0; j < 3; ++j) {
                if (b[k][j] == 0ULL) {
                    continue;
                }
                c[i][j] = static_cast<u64>((c[i][j] + static_cast<u128>(aik) * b[k][j]) % MOD);
            }
        }
    }
    return c;
}

Vec3 mat_vec_mul(const Mat3& a, const Vec3& v) {
    Vec3 out{};
    for (int i = 0; i < 3; ++i) {
        u128 sum = 0;
        for (int j = 0; j < 3; ++j) {
            if (a[i][j] == 0ULL || v[j] == 0ULL) {
                continue;
            }
            sum += static_cast<u128>(a[i][j]) * v[j];
        }
        out[i] = static_cast<u64>(sum % MOD);
    }
    return out;
}

Mat3 mat_pow(Mat3 base, u64 exp) {
    Mat3 result{};
    for (int i = 0; i < 3; ++i) {
        result[i][i] = 1ULL;
    }
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = mat_mul(base, result);
        }
        e >>= 1ULL;
        if (e > 0ULL) {
            base = mat_mul(base, base);
        }
    }
    return result;
}

u64 F_value(u64 n, u64 x_raw) {
    if (n == 0ULL) {
        return 0ULL;
    }

    const u64 x = x_raw % MOD;
    const u64 x2 = static_cast<u64>((static_cast<u128>(x) * x) % MOD);

    // a_n = f_n * x^n obeys: a_{n+1} = x*a_n + x^2*a_{n-1}.
    // State v_n = [a_n, a_{n-1}, F_n]^T with transition:
    // [a_{n+1}]   [x  x^2 0] [a_n    ]
    // [a_n    ] = [1  0   0] [a_{n-1}]
    // [F_{n+1}]   [x  x^2 1] [F_n    ].
    const Mat3 t{{
        {{x, x2, 0ULL}},
        {{1ULL, 0ULL, 0ULL}},
        {{x, x2, 1ULL}},
    }};

    const Vec3 v1{{x, 0ULL, x}};  // n=1
    const Mat3 p = mat_pow(t, n - 1ULL);
    const Vec3 vn = mat_vec_mul(p, v1);
    return vn[2];
}

u64 brute_small(int n, int x) {
    std::vector<u64> fib(static_cast<std::size_t>(n + 1), 0ULL);
    if (n >= 1) {
        fib[1] = 1ULL;
    }
    for (int i = 2; i <= n; ++i) {
        fib[static_cast<std::size_t>(i)] = fib[static_cast<std::size_t>(i - 1)] + fib[static_cast<std::size_t>(i - 2)];
    }

    u64 p = 1ULL;
    u64 sum = 0ULL;
    for (int i = 0; i <= n; ++i) {
        sum = (sum + static_cast<u128>(fib[static_cast<std::size_t>(i)] % MOD) * p) % MOD;
        p = static_cast<u64>((static_cast<u128>(p) * static_cast<u64>(x)) % MOD);
    }
    return sum;
}

u64 solve_sum(u64 n, int x_max) {
    u64 ans = 0ULL;
    for (int x = 0; x <= x_max; ++x) {
        ans += F_value(n, static_cast<u64>(x));
        ans %= MOD;
    }
    return ans;
}

bool run_checkpoints() {
    if (F_value(7ULL, 11ULL) != 268357683ULL) {
        std::cerr << "Checkpoint failed: F_7(11)\n";
        return false;
    }
    if (F_value(20ULL, 7ULL) != brute_small(20, 7)) {
        std::cerr << "Checkpoint failed: matrix vs brute for n=20,x=7\n";
        return false;
    }
    if (solve_sum(20ULL, 10) != 1044074802100ULL) {
        std::cerr << "Checkpoint failed: sum_{x=0..10} F_20(x)\n";
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

    std::cout << solve_sum(options.n, options.x_max) << '\n';
    return 0;
}
