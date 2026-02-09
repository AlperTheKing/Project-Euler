#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

long double integrand_center_speed(const long double a, const long double b, const long double t) {
    const long double s = std::sin(t);
    const long double c = std::cos(t);
    const long double num = a * b * std::sqrt(a * a * c * c + b * b * s * s);
    const long double den = a * a * s * s + b * b * c * c;
    return num / den;
}

long double simpson(const long double fa, const long double fm, const long double fb, const long double a,
                    const long double b) {
    return (b - a) * (fa + 4 * fm + fb) / 6;
}

long double adaptive_simpson(long double (*f)(long double, long double, long double), const long double a_param,
                             const long double b_param, const long double a, const long double b,
                             const long double eps, const long double whole, const long double fa,
                             const long double fm, const long double fb, int depth) {
    const long double m = (a + b) / 2;
    const long double l = (a + m) / 2;
    const long double r = (m + b) / 2;

    const long double fl = f(a_param, b_param, l);
    const long double fr = f(a_param, b_param, r);

    const long double left = simpson(fa, fl, fm, a, m);
    const long double right = simpson(fm, fr, fb, m, b);
    const long double delta = left + right - whole;

    if (depth <= 0 || std::fabsl(delta) <= 15 * eps) {
        // Richardson extrapolation.
        return left + right + delta / 15;
    }
    return adaptive_simpson(f, a_param, b_param, a, m, eps / 2, left, fa, fl, fm, depth - 1) +
           adaptive_simpson(f, a_param, b_param, m, b, eps / 2, right, fm, fr, fb, depth - 1);
}

long double integrate_center_curve_length(const long double a, const long double b) {
    // From rolling kinematics: if r(t)=(a cos t, b sin t) is the contact point relative to the center,
    // and φ(t) is the tangent angle, then the center velocity is |C'(t)| = |φ'(t)| * |r(t)|.
    //
    // With r'(t)=(-a sin t, b cos t), we get:
    //   φ'(t) = (x y' - y x') / |r'(t)|^2 = ab / (a^2 sin^2 t + b^2 cos^2 t),
    //   |r(t)| = sqrt(a^2 cos^2 t + b^2 sin^2 t).
    //
    // Hence the arc-length integrand is:
    //   ab * sqrt(a^2 cos^2 t + b^2 sin^2 t) / (a^2 sin^2 t + b^2 cos^2 t).
    //
    // The integrand is π-periodic and symmetric, so one full turn is:
    //   C(a,b) = 4 * ∫[0, π/2] integrand(t) dt.
    const long double pi = acosl(-1.0L);
    const long double A = 0.0L;
    const long double B = pi / 2;
    const long double fa = integrand_center_speed(a, b, A);
    const long double fb = integrand_center_speed(a, b, B);
    const long double m = (A + B) / 2;
    const long double fm = integrand_center_speed(a, b, m);
    const long double whole = simpson(fa, fm, fb, A, B);
    const long double eps = 1e-13L;
    const long double quarter = adaptive_simpson(&integrand_center_speed, a, b, A, B, eps, whole, fa, fm, fb, 30);
    return 4 * quarter;
}

bool run_checkpoints() {
    const long double v = integrate_center_curve_length(2.0L, 4.0L);
    const long double expected = 21.38816906L;
    if (std::fabsl(v - expected) > 5e-9L) {
        std::cerr << std::setprecision(15) << "Checkpoint failed: C(2,4) got " << v << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    const long double c14 = integrate_center_curve_length(1.0L, 4.0L);
    const long double c34 = integrate_center_curve_length(3.0L, 4.0L);
    const long double ans = c14 + c34;
    std::cout << std::fixed << std::setprecision(8) << ans << '\n';
    return 0;
}

