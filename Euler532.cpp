#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

long double simpson(long double fa, long double fm, long double fb, long double a, long double b) {
    return (b - a) * (fa + 4 * fm + fb) / 6;
}

long double adaptive_simpson(long double (*f)(long double, long double), long double p, long double a, long double b,
                             long double eps, long double whole, long double fa, long double fm, long double fb,
                             int depth) {
    const long double m = (a + b) / 2;
    const long double l = (a + m) / 2;
    const long double r = (m + b) / 2;
    const long double fl = f(l, p);
    const long double fr = f(r, p);

    const long double left = simpson(fa, fl, fm, a, m);
    const long double right = simpson(fm, fr, fb, m, b);
    const long double delta = left + right - whole;
    if (depth <= 0 || std::fabsl(delta) <= 15 * eps) {
        return left + right + delta / 15;
    }
    return adaptive_simpson(f, p, a, m, eps / 2, left, fa, fl, fm, depth - 1) +
           adaptive_simpson(f, p, m, b, eps / 2, right, fm, fr, fb, depth - 1);
}

// Integrand after substituting u = sin(theta), where theta is colatitude:
// L(n) = (1/s) * ∫_0^{u0} sqrt(1 - (s*u)^2) / (1 - u^2) du,  s = sin(pi/n).
long double integrand_u(long double u, long double s) {
    const long double uu = u * u;
    const long double su = s * u;
    return std::sqrt(1.0L - su * su) / (1.0L - uu);
}

long double length_per_bot(int n) {
    const long double pi = acosl(-1.0L);
    const long double s = std::sin(pi / static_cast<long double>(n));
    const long double u0 = 0.999L;

    const long double a = 0.0L;
    const long double b = u0;
    const long double m = (a + b) / 2;
    const long double fa = integrand_u(a, s);
    const long double fb = integrand_u(b, s);
    const long double fm = integrand_u(m, s);
    const long double whole = simpson(fa, fm, fb, a, b);
    const long double eps = 1e-13L;
    const long double I = adaptive_simpson(&integrand_u, s, a, b, eps, whole, fa, fm, fb, 30);
    return I / s;
}

bool run_checkpoints() {
    const long double L3 = length_per_bot(3);
    // Given: per-bot length rounds to 2.84 for n=3.
    if (std::fabsl(L3 - 2.84L) > 0.01L) {
        std::cerr << std::setprecision(15) << "Checkpoint failed: L(3) got " << L3 << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    int lo = 3;
    int hi = 3;
    while (length_per_bot(hi) <= 1000.0L) {
        hi *= 2;
    }
    while (lo + 1 < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (length_per_bot(mid) > 1000.0L) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const int n = hi;
    const long double per = length_per_bot(n);
    const long double per_prev = length_per_bot(n - 1);
    if (!(per > 1000.0L && per_prev <= 1000.0L)) {
        std::cerr << std::setprecision(15) << "Threshold check failed: n=" << n << " L(n)=" << per
                  << " L(n-1)=" << per_prev << '\n';
        return 2;
    }
    const long double total = static_cast<long double>(n) * per;
    std::cout << std::fixed << std::setprecision(2) << total << '\n';
    return 0;
}
