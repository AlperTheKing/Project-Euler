#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

using ld = long double;
using cd = std::complex<ld>;

static constexpr ld PI = 3.141592653589793238462643383279502884L;

static ld abs_one_minus_q_pow(ld s, ld t, int n) {
    ld rn = std::pow(s, static_cast<ld>(n));
    ld ang = static_cast<ld>(n) * t;
    cd qn = std::polar(rn, ang);
    return std::abs(cd(1, 0) - qn);
}

static std::pair<ld, ld> equations(ld s, ld t) {
    ld a1 = abs_one_minus_q_pow(s, t, 1);
    ld a7 = abs_one_minus_q_pow(s, t, 7);
    ld a8 = abs_one_minus_q_pow(s, t, 8);

    ld f7 = (1 + s) * a7 - (1 + std::pow(s, 7)) * a1;
    ld f8 = (1 + s) * a8 - (1 + std::pow(s, 8)) * a1;
    return {f7, f8};
}

static ld circular_triangle_area(ld r1, ld r2, ld r3) {
    ld a = r2 + r3;
    ld b = r1 + r3;
    ld c = r1 + r2;

    ld p = (a + b + c) / 2;
    ld tri = std::sqrt(std::max<ld>(0, p * (p - a) * (p - b) * (p - c)));

    ld A1 = std::acos((b * b + c * c - a * a) / (2 * b * c));
    ld A2 = std::acos((a * a + c * c - b * b) / (2 * a * c));
    ld A3 = std::acos((a * a + b * b - c * c) / (2 * a * b));

    return tri - 0.5L * (r1 * r1 * A1 + r2 * r2 * A2 + r3 * r3 * A3);
}

int main() {
    ld s = 0.9L;
    ld t = 0.83L;

    for (int it = 0; it < 80; ++it) {
        auto [f1, f2] = equations(s, t);

        ld h = 1e-12L;
        auto [f1s, f2s] = equations(s + h, t);
        auto [f1t, f2t] = equations(s, t + h);

        ld j11 = (f1s - f1) / h;
        ld j21 = (f2s - f2) / h;
        ld j12 = (f1t - f1) / h;
        ld j22 = (f2t - f2) / h;

        ld det = j11 * j22 - j12 * j21;
        assert(std::fabsl(det) > 1e-24L);

        ld ds = (-f1 * j22 + f2 * j12) / det;
        ld dt = (-j11 * f2 + j21 * f1) / det;

        s += ds;
        t += dt;

        while (t < 0) t += 2 * PI;
        while (t >= 2 * PI) t -= 2 * PI;

        if (std::fabsl(ds) + std::fabsl(dt) < 1e-20L) break;
    }

    auto [rf1, rf2] = equations(s, t);
    assert(std::fabsl(rf1) < 1e-12L);
    assert(std::fabsl(rf2) < 1e-12L);

    ld d = (1 + s) / abs_one_minus_q_pow(s, t, 1);
    cd q = std::polar(s, t);

    for (int i = 0; i < 20; ++i) {
        for (int j = i + 1; j < 40; ++j) {
            ld ri = std::pow(s, static_cast<ld>(i));
            ld rj = std::pow(s, static_cast<ld>(j));
            cd zi = d * std::pow(q, i);
            cd zj = d * std::pow(q, j);
            ld gap = std::abs(zi - zj) - (ri + rj);

            int diff = j - i;
            if (diff == 1 || diff == 7 || diff == 8) {
                assert(std::fabsl(gap) < 1e-10L);
            } else {
                assert(gap > -1e-10L);
            }
        }
    }

    ld A = circular_triangle_area(1, s, std::pow(s, 8));
    ld B = circular_triangle_area(1, std::pow(s, 7), std::pow(s, 8));
    ld total = (A + B) / (1 - s * s);

    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(total) << '\n';
    return 0;
}
