#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

constexpr std::uint64_t kMod = 1'123'581'313ULL;

using Mat2 = std::array<std::array<std::uint64_t, 2>, 2>;
using Vec2 = std::array<std::uint64_t, 2>;

Mat2 mat_mul(const Mat2& a, const Mat2& b) {
    Mat2 c{{{0, 0}, {0, 0}}};
    for (int i = 0; i < 2; ++i) {
        for (int k = 0; k < 2; ++k) {
            if (a[i][k] == 0) {
                continue;
            }
            for (int j = 0; j < 2; ++j) {
                c[i][j] = (c[i][j] + a[i][k] * b[k][j]) % kMod;
            }
        }
    }
    return c;
}

Mat2 mat_pow(Mat2 base, std::uint64_t exp) {
    Mat2 res{{{1, 0}, {0, 1}}};
    while (exp > 0) {
        if (exp & 1ULL) {
            res = mat_mul(res, base);
        }
        base = mat_mul(base, base);
        exp >>= 1ULL;
    }
    return res;
}

Vec2 mat_vec_mul(const Mat2& a, const Vec2& v) {
    Vec2 r{0, 0};
    for (int i = 0; i < 2; ++i) {
        r[i] = (a[i][0] * v[0] + a[i][1] * v[1]) % kMod;
    }
    return r;
}

std::uint64_t A_value(std::uint64_t m, std::uint64_t n) {
    const Mat2 C{{{0, 1}, {3, 1}}};
    const Mat2 M{{{1, 1}, {3, 2}}};
    const Vec2 v0{0, 1};

    const Mat2 mm = mat_pow(M, m);
    const Mat2 cn = mat_pow(C, n);
    const Vec2 vn = mat_vec_mul(cn, v0);
    return (mm[0][0] * vn[0] + mm[0][1] * vn[1]) % kMod;
}

std::uint64_t solve(int k) {
    std::vector<std::uint64_t> fib(static_cast<std::size_t>(k) + 1, 0);
    fib[0] = 0;
    fib[1] = 1;
    for (int i = 2; i <= k; ++i) {
        fib[i] = fib[i - 1] + fib[i - 2];
    }

    const Mat2 C{{{0, 1}, {3, 1}}};
    const Mat2 M{{{1, 1}, {3, 2}}};
    const Vec2 v0{0, 1};

    std::vector<Vec2> col(static_cast<std::size_t>(k) + 1);
    std::vector<Mat2> row(static_cast<std::size_t>(k) + 1);

    for (int i = 2; i <= k; ++i) {
        row[i] = mat_pow(M, fib[i]);
        col[i] = mat_vec_mul(mat_pow(C, fib[i]), v0);
    }

    std::uint64_t ans = 0;
    for (int i = 2; i <= k; ++i) {
        for (int j = 2; j <= k; ++j) {
            const std::uint64_t val =
                (row[i][0][0] * col[j][0] + row[i][0][1] * col[j][1]) % kMod;
            ans += val;
            if (ans >= kMod) {
                ans -= kMod;
            }
        }
    }
    return ans;
}

void run_validations() {
    assert(A_value(1, 1) == 2ULL);
    assert(A_value(1, 2) == 5ULL);
    assert(A_value(2, 1) == 7ULL);
    assert(A_value(2, 2) == 16ULL);
    assert(solve(3) == 30ULL);
    assert(solve(5) == 10'396ULL);
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve(50) << '\n';
    return 0;
}
