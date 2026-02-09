#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = long long;

constexpr u64 kMod = 1000000000ULL;  // last 9 digits
constexpr int kDim = 11;

struct Options {
    i64 k_max = 1234567890123LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoll(tail);
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
        if (parse_i64_after_prefix(arg, "--k-max=", options.k_max)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.k_max >= 2;
}

i64 gcd4(i64 a, i64 b, i64 c, i64 d) {
    auto g2 = [](i64 x, i64 y) {
        while (y != 0) {
            const i64 t = x % y;
            x = y;
            y = t;
        }
        return x < 0 ? -x : x;
    };
    return g2(g2(a, b), g2(c, d));
}

struct Coefs {
    u64 a3 = 0;  // coefficient of q^3
    std::array<u64, 24> a2{};  // coefficient of q^2 by remainder
    std::array<u64, 24> a1{};  // coefficient of q   by remainder
    std::array<u64, 24> a0{};  // constant by remainder
};

Coefs build_coefficients() {
    Coefs out;

    for (int ra = 1; ra <= 24; ++ra) {
        for (int rb = 1; rb <= 24; ++rb) {
            for (int rc = 1; rc <= 24; ++rc) {
                const i64 m = gcd4(
                    24,
                    36 + 6 * ra,
                    14 + 6 * ra + 2 * rb,
                    1 + ra + rb + rc);

                out.a3 = (out.a3 + static_cast<u64>(m)) % kMod;
                for (int rem = 0; rem < 24; ++rem) {
                    const int ea = (ra <= rem ? 1 : 0);
                    const int eb = (rb <= rem ? 1 : 0);
                    const int ec = (rc <= rem ? 1 : 0);
                    const u64 u = static_cast<u64>(m) % kMod;
                    out.a2[static_cast<std::size_t>(rem)] =
                        (out.a2[static_cast<std::size_t>(rem)] + u * static_cast<u64>(ea + eb + ec)) % kMod;
                    out.a1[static_cast<std::size_t>(rem)] =
                        (out.a1[static_cast<std::size_t>(rem)] + u * static_cast<u64>(ea * eb + ea * ec + eb * ec)) % kMod;
                    out.a0[static_cast<std::size_t>(rem)] =
                        (out.a0[static_cast<std::size_t>(rem)] + u * static_cast<u64>(ea * eb * ec)) % kMod;
                }
            }
        }
    }

    return out;
}

u64 s_mod(i64 n, const Coefs& coef) {
    const i64 q = n / 24;
    const int rem = static_cast<int>(n % 24);
    const u64 qq = static_cast<u64>(q % static_cast<i64>(kMod));
    const u64 q2 = static_cast<u64>((__uint128_t)qq * qq % kMod);
    const u64 q3 = static_cast<u64>((__uint128_t)q2 * qq % kMod);

    u64 ans = 0;
    ans = (ans + static_cast<u64>((__uint128_t)coef.a3 * q3 % kMod)) % kMod;
    ans = (ans + static_cast<u64>((__uint128_t)coef.a2[static_cast<std::size_t>(rem)] * q2 % kMod)) % kMod;
    ans = (ans + static_cast<u64>((__uint128_t)coef.a1[static_cast<std::size_t>(rem)] * qq % kMod)) % kMod;
    ans = (ans + coef.a0[static_cast<std::size_t>(rem)]) % kMod;
    return ans;
}

struct Mat {
    std::array<std::array<u64, kDim>, kDim> a{};
};

Mat mat_identity() {
    Mat m;
    for (int i = 0; i < kDim; ++i) {
        m.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1ULL;
    }
    return m;
}

Mat mat_mul(const Mat& x, const Mat& y) {
    Mat z;
    for (int i = 0; i < kDim; ++i) {
        for (int k = 0; k < kDim; ++k) {
            const u64 xik = x.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            if (xik == 0ULL) {
                continue;
            }
            for (int j = 0; j < kDim; ++j) {
                const u64 ykj = y.a[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
                if (ykj == 0ULL) {
                    continue;
                }
                z.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                    (z.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                     static_cast<u64>((__uint128_t)xik * ykj % kMod)) %
                    kMod;
            }
        }
    }
    return z;
}

Mat mat_pow(Mat base, i64 exp) {
    Mat res = mat_identity();
    i64 e = exp;
    while (e > 0) {
        if (e & 1LL) {
            res = mat_mul(base, res);
        }
        base = mat_mul(base, base);
        e >>= 1LL;
    }
    return res;
}

std::array<u64, kDim> mat_vec_mul(const Mat& m, const std::array<u64, kDim>& v) {
    std::array<u64, kDim> out{};
    for (int i = 0; i < kDim; ++i) {
        u64 s = 0ULL;
        for (int j = 0; j < kDim; ++j) {
            const u64 mij = m.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            if (mij == 0ULL) {
                continue;
            }
            s = (s + static_cast<u64>((__uint128_t)mij * v[static_cast<std::size_t>(j)] % kMod)) % kMod;
        }
        out[static_cast<std::size_t>(i)] = s;
    }
    return out;
}

Mat step_matrix(int carry, int rem, const Coefs& coef) {
    Mat m{};
    const u64 c = static_cast<u64>(carry) % kMod;
    const u64 c2 = static_cast<u64>((__uint128_t)c * c % kMod);
    const u64 c3 = static_cast<u64>((__uint128_t)c2 * c % kMod);

    // Basis:
    // 0:1, 1:q, 2:q1, 3:q^2, 4:q q1, 5:q1^2, 6:q^3, 7:q^2 q1, 8:q q1^2, 9:q1^3, 10:sum
    m.a[0][0] = 1;

    m.a[1][2] = 1;                // q' = q1
    m.a[2][0] = c;                // q1' = q + q1 + c
    m.a[2][1] = 1;
    m.a[2][2] = 1;

    m.a[3][5] = 1;                // q'^2 = q1^2

    m.a[4][2] = c;                // q' q1' = q q1 + q1^2 + c q1
    m.a[4][4] = 1;
    m.a[4][5] = 1;

    m.a[5][0] = c2;               // q1'^2
    m.a[5][1] = (2 * c) % kMod;
    m.a[5][2] = (2 * c) % kMod;
    m.a[5][3] = 1;
    m.a[5][4] = 2;
    m.a[5][5] = 1;

    m.a[6][9] = 1;                // q'^3 = q1^3

    m.a[7][5] = c;                // q'^2 q1'
    m.a[7][8] = 1;
    m.a[7][9] = 1;

    m.a[8][2] = c2;               // q' q1'^2
    m.a[8][4] = (2 * c) % kMod;
    m.a[8][5] = (2 * c) % kMod;
    m.a[8][7] = 1;
    m.a[8][8] = 2;
    m.a[8][9] = 1;

    m.a[9][0] = c3;               // q1'^3
    m.a[9][1] = (3 * c2) % kMod;
    m.a[9][2] = (3 * c2) % kMod;
    m.a[9][3] = (3 * c) % kMod;
    m.a[9][4] = (6 * c) % kMod;
    m.a[9][5] = (3 * c) % kMod;
    m.a[9][6] = 1;
    m.a[9][7] = 3;
    m.a[9][8] = 3;
    m.a[9][9] = 1;

    // sum' = sum + S(F_k), and S(F_k)=a3*q^3 + a2(rem)*q^2 + a1(rem)*q + a0(rem)
    m.a[10][0] = coef.a0[static_cast<std::size_t>(rem)];
    m.a[10][1] = coef.a1[static_cast<std::size_t>(rem)];
    m.a[10][3] = coef.a2[static_cast<std::size_t>(rem)];
    m.a[10][6] = coef.a3;
    m.a[10][10] = 1;

    return m;
}

u64 solve(const i64 k_max, const Coefs& coef) {
    if (k_max < 2) {
        return 0;
    }

    std::array<int, 24> fib_mod24{};
    fib_mod24[0] = 0;
    fib_mod24[1] = 1;
    for (int i = 2; i < 24; ++i) {
        fib_mod24[static_cast<std::size_t>(i)] =
            (fib_mod24[static_cast<std::size_t>(i - 1)] + fib_mod24[static_cast<std::size_t>(i - 2)]) % 24;
    }

    std::array<int, 24> carry{};
    for (int i = 0; i < 24; ++i) {
        carry[static_cast<std::size_t>(i)] =
            (fib_mod24[static_cast<std::size_t>(i)] + fib_mod24[static_cast<std::size_t>((i + 1) % 24)]) / 24;
    }

    std::array<Mat, 24> step{};
    for (int phase = 0; phase < 24; ++phase) {
        step[static_cast<std::size_t>(phase)] =
            step_matrix(carry[static_cast<std::size_t>(phase)],
                        fib_mod24[static_cast<std::size_t>(phase)],
                        coef);
    }

    // We apply steps for k = 2..k_max, so number of steps is k_max-1.
    const i64 steps = k_max - 1;
    const int start_phase = 2 % 24;

    Mat block = mat_identity();
    for (int t = 0; t < 24; ++t) {
        const int phase = (start_phase + t) % 24;
        block = mat_mul(step[static_cast<std::size_t>(phase)], block);
    }

    std::array<u64, kDim> state{};
    // F_2 = 1 = 24*0 + 1, F_3 = 2 = 24*0 + 2 => q_2=0, q_3=0.
    state[0] = 1;  // constant term
    state[1] = 0;
    state[2] = 0;
    state[3] = 0;
    state[4] = 0;
    state[5] = 0;
    state[6] = 0;
    state[7] = 0;
    state[8] = 0;
    state[9] = 0;
    state[10] = 0;  // accumulated sum

    const i64 full_blocks = steps / 24;
    const int rem_steps = static_cast<int>(steps % 24);

    if (full_blocks > 0) {
        const Mat block_pow = mat_pow(block, full_blocks);
        state = mat_vec_mul(block_pow, state);
    }

    for (int t = 0; t < rem_steps; ++t) {
        const int phase = (start_phase + t) % 24;
        state = mat_vec_mul(step[static_cast<std::size_t>(phase)], state);
    }

    return state[10] % kMod;
}

bool run_checkpoints(const Coefs& coef) {
    if (s_mod(10, coef) != 1972ULL) {
        std::cerr << "Checkpoint failed: S(10)\n";
        return false;
    }
    if (s_mod(10000, coef) != (2024258331114ULL % kMod)) {
        std::cerr << "Checkpoint failed: S(10000) mod 1e9\n";
        return false;
    }
    // Exact S(10000) is given; recompute directly in 128-bit for full check.
    {
        const i64 n = 10000;
        const i64 q = n / 24;
        const int rem = static_cast<int>(n % 24);
        __int128 total = 0;
        for (int a = 1; a <= 24; ++a) {
            const i64 ca = q + (a <= rem ? 1 : 0);
            for (int b = 1; b <= 24; ++b) {
                const i64 cb = q + (b <= rem ? 1 : 0);
                for (int c = 1; c <= 24; ++c) {
                    const i64 cc = q + (c <= rem ? 1 : 0);
                    const i64 m = gcd4(
                        24,
                        36 + 6 * a,
                        14 + 6 * a + 2 * b,
                        1 + a + b + c);
                    total += static_cast<__int128>(ca) * cb * cc * m;
                }
            }
        }
        const i64 exact = static_cast<i64>(total);
        if (exact != 2024258331114LL) {
            std::cerr << "Checkpoint failed: exact S(10000)\n";
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const Coefs coef = build_coefficients();
    if (options.run_checkpoints && !run_checkpoints(coef)) {
        return 2;
    }

    const u64 answer = solve(options.k_max, coef);
    std::cout << std::setfill('0') << std::setw(9) << answer << '\n';
    return 0;
}
