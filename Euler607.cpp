#include <cmath>
#include <iomanip>
#include <iostream>

// Project Euler 607: layered medium with parallel boundaries; use Snell's law (sin(theta)/v constant)
// after rotating coordinates so the marsh boundaries become horizontal.

int main() {
    constexpr long double AtoB = 100.0L;
    const long double rt2 = std::sqrt(2.0L);

    // In rotated coords: v = (x-y)/sqrt(2) is the normal direction to the marsh boundaries.
    const long double d0 = 25.0L * (rt2 - 1.0L); // outside thickness in v on each side
    const long double dv[7] = {d0, 10.0L, 10.0L, 10.0L, 10.0L, 10.0L, d0};
    const long double sp[7] = {10.0L, 9.0L, 8.0L, 7.0L, 6.0L, 5.0L, 10.0L};

    const long double U = AtoB / rt2; // required delta-u in rotated coords

    auto horiz = [&](long double k) {
        long double s = 0.0L;
        for (int i = 0; i < 7; ++i) {
            const long double t = k * sp[i];
            const long double c = std::sqrt(1.0L - t * t);
            s += dv[i] * t / c;
        }
        return s;
    };

    auto time_needed = [&](long double k) {
        long double s = 0.0L;
        for (int i = 0; i < 7; ++i) {
            const long double t = k * sp[i];
            const long double c = std::sqrt(1.0L - t * t);
            s += dv[i] / (sp[i] * c);
        }
        return s;
    };

    // Binary search k in (0, 1/max_speed).
    long double lo = 0.0L, hi = 0.1L;
    for (int it = 0; it < 200; ++it) {
        const long double mid = (lo + hi) / 2;
        if (horiz(mid) < U) lo = mid;
        else hi = mid;
    }

    // Statement sanity: direct-east path time.
    const long double direct = (dv[0] * rt2) / 10.0L +
                               (10.0L * rt2) / 9.0L + (10.0L * rt2) / 8.0L + (10.0L * rt2) / 7.0L +
                               (10.0L * rt2) / 6.0L + (10.0L * rt2) / 5.0L +
                               (dv[6] * rt2) / 10.0L;
    if (std::fabsl(direct - 13.4738L) > 1e-3L) {
        std::cerr << "Sanity check failed\n";
        return 1;
    }

    const long double ans = time_needed((lo + hi) / 2);
    std::cout << std::fixed << std::setprecision(10) << (double)ans << "\n";
    return 0;
}

