#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {

long double pi_ld() {
    return acosl(-1.0L);
}

long double clamp_ld(long double x, long double lo, long double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

long double lambda_regular(int n) {
    long double theta = pi_ld() / static_cast<long double>(n);
    long double tan_theta = tanl(theta);
    int K = -1;
    for (int k = 0; k <= n; ++k) {
        long double ktheta = static_cast<long double>(k) * theta;
        long double val = sinl(ktheta) - (static_cast<long double>(k + n) * tan_theta * cosl(ktheta));
        if (val < 0.0L) {
            K = k;
        } else {
            break;
        }
    }
    if (K < 0) K = 0;

    long double ktheta = static_cast<long double>(K) * theta;
    long double denom = static_cast<long double>(K + n) * tan_theta;
    long double arg = 2.0L * sinl(ktheta) / denom - cosl(ktheta);
    arg = clamp_ld(arg, -1.0L, 1.0L);

    long double alpha = 0.5L * (ktheta + acosl(arg));
    return 1.0L / cosl(alpha);
}

long double lambda_circle() {
    long double lo = 0.0L;
    long double hi = pi_ld() / 2.0L - 1e-12L;
    auto f = [](long double x) { return x + pi_ld() - tanl(x); };
    long double flo = f(lo);
    long double fhi = f(hi);
    if (!(flo > 0.0L && fhi < 0.0L)) {
        return std::numeric_limits<long double>::quiet_NaN();
    }
    for (int iter = 0; iter < 200; ++iter) {
        long double mid = 0.5L * (lo + hi);
        long double fm = f(mid);
        if (fm > 0.0L) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    long double mu = 0.5L * (lo + hi);
    return 1.0L / cosl(mu);
}

}  // namespace

int main() {
    const long double expected_circle = 4.60333885L;
    const long double expected_square = 5.78859314L;
    const long double tol = 5e-8L;

    long double circle = lambda_circle();
    if (!std::isfinite(circle) || fabsl(circle - expected_circle) > tol) {
        std::cerr << "Validation failed for circle: got "
                  << std::fixed << std::setprecision(12) << circle
                  << " expected " << expected_circle << "\n";
        return 1;
    }

    long double square = lambda_regular(4);
    if (fabsl(square - expected_square) > tol) {
        std::cerr << "Validation failed for square: got "
                  << std::fixed << std::setprecision(12) << square
                  << " expected " << expected_square << "\n";
        return 1;
    }

    long double hexagon = lambda_regular(6);
    std::cout.setf(std::ios::fixed);
    std::cout.precision(8);
    std::cout << static_cast<double>(hexagon) << "\n";
    return 0;
}
