#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

long double curve_area(const long double v) {
    // Integral of x(t) * y'(t) over t in [0,1].
    return 0.5L + (3.0L * v) / 5.0L - (3.0L * v * v) / 20.0L;
}

long double curve_speed(const long double t, const long double v) {
    const long double dx = 6.0L * (v - 1.0L) * t + 3.0L * (2.0L - 3.0L * v) * t * t;
    const long double dy = 3.0L * v + 2.0L * (3.0L - 6.0L * v) * t + 3.0L * (3.0L * v - 2.0L) * t * t;
    return std::sqrt(dx * dx + dy * dy);
}

long double adaptive_simpson(const long double v,
                             const long double a,
                             const long double b,
                             const long double eps,
                             const long double whole,
                             const long double fa,
                             const long double fb,
                             const long double fm,
                             const int depth) {
    const long double m = (a + b) * 0.5L;
    const long double lm = (a + m) * 0.5L;
    const long double rm = (m + b) * 0.5L;

    const long double flm = curve_speed(lm, v);
    const long double frm = curve_speed(rm, v);

    const long double left = (m - a) * (fa + 4.0L * flm + fm) / 6.0L;
    const long double right = (b - m) * (fm + 4.0L * frm + fb) / 6.0L;
    const long double combined = left + right;

    if (depth <= 0 || std::fabsl(combined - whole) <= 15.0L * eps) {
        return combined + (combined - whole) / 15.0L;
    }

    return adaptive_simpson(v, a, m, eps * 0.5L, left, fa, fm, flm, depth - 1) +
           adaptive_simpson(v, m, b, eps * 0.5L, right, fm, fb, frm, depth - 1);
}

long double curve_length(const long double v) {
    const long double a = 0.0L;
    const long double b = 1.0L;
    const long double fa = curve_speed(a, v);
    const long double fb = curve_speed(b, v);
    const long double m = 0.5L;
    const long double fm = curve_speed(m, v);
    const long double whole = (b - a) * (fa + 4.0L * fm + fb) / 6.0L;
    return adaptive_simpson(v, a, b, 1e-18L, whole, fa, fb, fm, 30);
}

long double solve_v() {
    // From area equation:
    //   1/2 + 3v/5 - 3v^2/20 = pi/4
    // => v = 2 +/- sqrt(66 - 15*pi) / 3
    const long double disc = 66.0L - 15.0L * std::acos(-1.0L);
    return 2.0L - std::sqrt(disc) / 3.0L;
}

bool run_checkpoints() {
    if (std::fabsl(curve_area(0.0L) - 0.5L) > 1e-18L) {
        std::cerr << "Checkpoint failed: area(v=0)\n";
        return false;
    }

    if (std::fabsl(curve_length(0.0L) - std::sqrt(2.0L)) > 1e-12L) {
        std::cerr << "Checkpoint failed: length(v=0)\n";
        return false;
    }

    const long double v = solve_v();
    if (!(v > 0.0L && v < 1.0L)) {
        std::cerr << "Checkpoint failed: v range\n";
        return false;
    }

    const long double target_area = std::acos(-1.0L) / 4.0L;
    if (std::fabsl(curve_area(v) - target_area) > 1e-18L) {
        std::cerr << "Checkpoint failed: area match\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const long double pi = std::acos(-1.0L);
    const long double v = solve_v();
    const long double L = curve_length(v);
    const long double percent = 100.0L * (L - pi / 2.0L) / (pi / 2.0L);

    std::cout << std::fixed << std::setprecision(10) << percent << '\n';
    return 0;
}
