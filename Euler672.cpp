#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1'117'117'717ULL;
constexpr int DIM = 5;

struct Mat {
    std::array<std::array<u64, DIM>, DIM> a{};

    static Mat identity() {
        Mat m;
        for (int i = 0; i < DIM; ++i) {
            m.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1ULL;
        }
        return m;
    }
};

Mat multiply(const Mat& x, const Mat& y) {
    Mat z;
    for (int i = 0; i < DIM; ++i) {
        for (int k = 0; k < DIM; ++k) {
            const u64 xik = x.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
            if (xik == 0ULL) {
                continue;
            }
            for (int j = 0; j < DIM; ++j) {
                const u64 ykj = y.a[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
                if (ykj == 0ULL) {
                    continue;
                }
                u64& cell = z.a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                cell = static_cast<u64>((static_cast<u128>(cell) +
                                         static_cast<u128>(xik) * static_cast<u128>(ykj)) %
                                        MOD);
            }
        }
    }
    return z;
}

Mat power(Mat base, u64 exp) {
    Mat result = Mat::identity();
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = multiply(base, result);
        }
        exp >>= 1ULL;
        if (exp > 0ULL) {
            base = multiply(base, base);
        }
    }
    return result;
}

std::array<u64, DIM> apply_matrix(const Mat& m, const std::array<u64, DIM>& v) {
    std::array<u64, DIM> out{};
    for (int i = 0; i < DIM; ++i) {
        u64 s = 0ULL;
        for (int j = 0; j < DIM; ++j) {
            s = static_cast<u64>((static_cast<u128>(s) +
                                  static_cast<u128>(m.a[static_cast<std::size_t>(i)]
                                                        [static_cast<std::size_t>(j)]) *
                                      static_cast<u128>(v[static_cast<std::size_t>(j)])) %
                                 MOD);
        }
        out[static_cast<std::size_t>(i)] = s;
    }
    return out;
}

u64 slow_g(u64 n) {
    u64 added = 0ULL;
    while (n > 1ULL) {
        if (n % 7ULL == 0ULL) {
            n /= 7ULL;
        } else {
            ++n;
            ++added;
        }
    }
    return added;
}

Mat digit_matrix(int d, const std::array<int, 7>& weight, const std::array<int, 7>& pref) {
    Mat m;

    m.a[0][0] = 1ULL;

    m.a[1][0] = 6ULL;
    m.a[1][1] = 7ULL;
    m.a[1][4] = static_cast<u64>(d);

    m.a[2][0] = 15ULL;
    m.a[2][1] = 21ULL;
    m.a[2][2] = 7ULL;
    m.a[2][3] = static_cast<u64>(d);
    m.a[2][4] = static_cast<u64>(pref[static_cast<std::size_t>(d)]);

    m.a[3][3] = 1ULL;
    m.a[3][4] = static_cast<u64>(weight[static_cast<std::size_t>(d)]);

    m.a[4][4] = 1ULL;

    return m;
}

Mat sequence_matrix(const std::string& seq, const std::array<int, 7>& weight,
                    const std::array<int, 7>& pref) {
    Mat all = Mat::identity();
    for (char ch : seq) {
        const int d = ch - '0';
        const Mat md = digit_matrix(d, weight, pref);
        all = multiply(md, all);
    }
    return all;
}

u64 solve_H(u64 k) {
    assert(k % 10ULL == 0ULL);
    assert(k >= 10ULL);

    const std::array<int, 7> weight{6, 5, 4, 3, 2, 1, 0};
    std::array<int, 7> pref{};
    int run = 0;
    for (int d = 0; d < 7; ++d) {
        pref[static_cast<std::size_t>(d)] = run;
        run += weight[static_cast<std::size_t>(d)];
    }

    const std::string blockA = "4311623550";
    const std::string blockB = "431162355";
    const std::string remSingle = "31162355";
    const std::string remFirst = "311623550";

    const Mat MA = sequence_matrix(blockA, weight, pref);
    const Mat MB = sequence_matrix(blockB, weight, pref);
    const Mat MSingle = sequence_matrix(remSingle, weight, pref);
    const Mat MFirst = sequence_matrix(remFirst, weight, pref);

    const u64 t = k / 10ULL;

    std::array<u64, DIM> v{};
    v[0] = 1ULL;
    v[1] = 3ULL;
    v[2] = 12ULL;
    v[3] = 2ULL;
    v[4] = 1ULL;

    if (t == 1ULL) {
        v = apply_matrix(MSingle, v);
    } else {
        v = apply_matrix(MFirst, v);
        if (t > 2ULL) {
            const Mat mid = power(MA, t - 2ULL);
            v = apply_matrix(mid, v);
        }
        v = apply_matrix(MB, v);
    }

    return v[2] % MOD;
}

}  // namespace

int main() {
    assert(slow_g(125ULL) == 8ULL);
    assert(slow_g(1'000ULL) == 9ULL);
    assert(slow_g(10'000ULL) == 21ULL);

    assert(solve_H(10ULL) == 690'409'338ULL);

    std::cout << solve_H(1'000'000'000ULL) << "\n";
    return 0;
}
