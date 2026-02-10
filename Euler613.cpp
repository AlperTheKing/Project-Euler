#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Project Euler 613: average the angular measure that sees the hypotenuse over all interior points.

using ld = long double;

static void gauss_legendre_01(int n, std::vector<ld>& x, std::vector<ld>& w) {
    x.assign(n, 0);
    w.assign(n, 0);
    const int m = (n + 1) / 2;
    const ld pi = acosl(-1.0L);
    const ld eps = 1e-18L;

    for (int i = 0; i < m; ++i) {
        ld z = cosl(pi * (i + 0.75L) / (n + 0.5L));
        for (int it = 0; it < 100; ++it) {
            ld p1 = 1.0L, p2 = 0.0L;
            for (int j = 1; j <= n; ++j) {
                const ld p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
            }
            const ld pp = n * (z * p1 - p2) / (z * z - 1.0L);
            const ld z1 = z;
            z = z1 - p1 / pp;
            if (z == z1 || fabsl(z - z1) < eps) break;
        }

        // Recompute derivative at converged root for the weight.
        ld p1 = 1.0L, p2 = 0.0L;
        for (int j = 1; j <= n; ++j) {
            const ld p3 = p2;
            p2 = p1;
            p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
        }
        const ld pp = n * (z * p1 - p2) / (z * z - 1.0L);
        const ld ww = 2.0L / ((1.0L - z * z) * pp * pp);

        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = ww;
        w[n - 1 - i] = ww;
    }

    for (int i = 0; i < n; ++i) {
        x[i] = 0.5L * (x[i] + 1.0L);
        w[i] *= 0.5L;
    }
}

static ld probability_hypotenuse(int n) {
    std::vector<ld> x, w;
    gauss_legendre_01(n, x, w);

    const ld pi = acosl(-1.0L);
    ld integral = 0.0L; // integral of angle over the triangle area

    for (int i = 0; i < n; ++i) {
        const ld u = x[i];
        const ld wu = w[i];
        const ld omu = 1.0L - u;

        const ld px = 40.0L * u;
        const ld jac = 1200.0L * omu; // dx dy = 1200 (1-u) du dv

        for (int j = 0; j < n; ++j) {
            const ld v = x[j];
            const ld wv = w[j];
            const ld py = 30.0L * omu * v;

            const ld ux = 40.0L - px, uy = -py;
            const ld vx = -px, vy = 30.0L - py;
            const ld cross = ux * vy - uy * vx;
            const ld dot = ux * vx + uy * vy;
            const ld theta = atan2l(fabsl(cross), dot);

            integral += wu * wv * theta * jac;
        }
    }

    // probability = (1 / (2*pi*Area)) * integral, and Area = 600, so denom = 1200*pi.
    return integral / (1200.0L * pi);
}

int main() {
    const ld p512 = probability_hypotenuse(512);
    const ld p1024 = probability_hypotenuse(1024);
    assert(fabsl(p1024 - p512) < 5e-12L);

    std::cout << std::fixed << std::setprecision(10) << p1024 << "\n";
    return 0;
}
