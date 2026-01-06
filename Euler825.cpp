#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using int64 = long long;
using i128 = __int128_t;
using ld = long double;

namespace {

ld digamma_ld(ld x) {
    if (!(x > 0.0L) || !std::isfinite(static_cast<double>(x))) {
        return std::numeric_limits<ld>::quiet_NaN();
    }

    ld res = 0.0L;
    while (x < 10.0L) {
        res -= 1.0L / x;
        x += 1.0L;
    }

    const ld inv = 1.0L / x;
    const ld inv2 = inv * inv;

    ld series = inv2 * (-1.0L / 12.0L
                 + inv2 * (1.0L / 120.0L
                 + inv2 * (-1.0L / 252.0L
                 + inv2 * (1.0L / 240.0L
                 + inv2 * (-1.0L / 132.0L
                 + inv2 * (691.0L / 32760.0L))))));

    res += std::logl(x) - 0.5L * inv + series;
    return res;
}

inline ld i128_to_ld(i128 v) {
    return static_cast<ld>(v);
}

} // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const unsigned long long N = 100000000000000ULL;

    const ld sqrt3 = std::sqrtl(3.0L);
    const ld c = (3.0L + sqrt3) / 6.0L;

    const int K = 60;
    std::vector<i128> num(K + 1, 0);
    std::vector<i128> den(K + 1, 0);

    num[2] = 7;
    num[3] = 35;
    num[4] = 121;

    den[2] = 11;
    den[3] = 73;
    den[4] = 395;
    den[5] = 1933;

    for (int n = 5; n <= K; ++n) {
        num[n] = 3 * num[n - 1] + 3 * num[n - 2] - num[n - 3];
    }
    for (int n = 6; n <= K; ++n) {
        den[n] = 8 * den[n - 1] - 18 * den[n - 2] + 8 * den[n - 3] - den[n - 4];
    }

    ld E = 0.0L;
    for (int n = 2; n <= K; ++n) {
        const ld s = i128_to_ld(num[n]) / i128_to_ld(den[n]);
        E += s - 1.0L / (static_cast<ld>(n) - c);
    }

#ifdef VALIDATE
    {
        const ld s2 = i128_to_ld(num[2]) / i128_to_ld(den[2]);
        const ld expected_s2 = 7.0L / 11.0L;
        if (std::fabsl(s2 - expected_s2) > 1e-18L) {
            throw std::runtime_error("Validation failed: S(2) mismatch");
        }

        ld t10 = 0.0L;
        for (int n = 2; n <= 10; ++n) {
            t10 += i128_to_ld(num[n]) / i128_to_ld(den[n]);
        }
        const ld expected_t10 = 2.38235282L;
        if (std::fabsl(t10 - expected_t10) > 5e-8L) {
            throw std::runtime_error("Validation failed: T(10) mismatch");
        }

        const int ncheck = 50;
        const ld s = i128_to_ld(num[ncheck]) / i128_to_ld(den[ncheck]);
        const ld delta = 1.0L / s - static_cast<ld>(ncheck);
        if (std::fabsl(delta + c) > 1e-12L) {
            throw std::runtime_error("Validation failed: asymptotic shift mismatch");
        }
    }
#endif

    const ld x_big = static_cast<ld>(N) + 1.0L - c;
    const ld x_small = 2.0L - c;

    const ld ans = digamma_ld(x_big) - digamma_ld(x_small) + E;

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(8) << ans << "\n";
    return 0;
}
