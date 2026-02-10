#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

// Project Euler 587: Concave Triangle
//
// Scale radius to 1. Then the bounding rectangle is [0, 2n] x [0, 2] and the leftmost
// circle is centered at (1,1) with radius 1. The L-section is the region in the unit
// square [0,1] x [0,1] outside that circle, so:
//   A_L = 1 - pi/4.
//
// The diagonal line from (0,0) to (2n,2) has equation y = x/n.
// The concave triangle is the part of the L-section lying below this line.
//
// Circle arc in the unit square is given by (x-1)^2 + (y-1)^2 = 1 with y <= 1:
//   y_arc(x) = 1 - sqrt(2x - x^2),  x in [0,1].
//
// The line intersects the arc at x = x0 where x0/n = y_arc(x0). Solving gives:
//   x0 = n (n + 1 - sqrt(2n)) / (n^2 + 1)   (the root in (0,1)).
//
// The area below the min(line, arc) is:
//   A_C(n) = integral_0^{x0} (x/n) dx + integral_{x0}^1 y_arc(x) dx.
// The second integral can be expressed via a quarter-circle segment.

using i64 = std::int64_t;

static long double pi_ld() {
    return acosl(-1.0L);
}

static long double concave_ratio(i64 n) {
    const long double nl = (long double)n;
    const long double s2n = sqrtl(2.0L * nl);
    const long double x0 = nl * (nl + 1.0L - s2n) / (nl * nl + 1.0L);

    const long double d = 1.0L - x0; // d in (0,1)
    const long double root = sqrtl(std::max(0.0L, 1.0L - d * d));

    // Integral pieces.
    const long double area_line = x0 * x0 / (2.0L * nl);
    const long double area_arc = d - 0.5L * (d * root + asinl(d));
    const long double area_concave = area_line + area_arc;

    const long double area_L = 1.0L - pi_ld() / 4.0L;
    return area_concave / area_L;
}

static i64 least_n(long double threshold) {
    for (i64 n = 1;; ++n) {
        if (concave_ratio(n) < threshold) return n;
    }
}

int main() {
    // Validation checks from the statement.
    {
        const long double r1 = concave_ratio(1);
        if (fabsl(r1 - 0.5L) > 1e-15L) {
            std::cerr << "Validation failed: ratio(1) got " << std::setprecision(18) << (double)r1 << "\n";
            return 1;
        }
    }
    {
        const long double r2 = concave_ratio(2);
        if (fabsl(r2 - 0.364626169291728L) > 1e-12L) {
            std::cerr << "Validation failed: ratio(2) got " << std::setprecision(18) << (double)r2 << "\n";
            return 1;
        }
    }
    {
        const i64 n10 = least_n(0.1L);
        if (n10 != 15) {
            std::cerr << "Validation failed: least_n(0.1) got " << n10 << "\n";
            return 1;
        }
    }

    std::cout << least_n(0.001L) << "\n";
    return 0;
}
