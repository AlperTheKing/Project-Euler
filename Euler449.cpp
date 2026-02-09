#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    long double a = 3.0L;
    long double b = 1.0L;
    bool run_checkpoints = true;
};

bool parse_ld_after_prefix(const std::string& arg, const std::string& prefix, long double& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stold(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_ld_after_prefix(arg, "--a=", options.a)) {
            continue;
        }
        if (parse_ld_after_prefix(arg, "--b=", options.b)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.a > 0.0L && options.b > 0.0L;
}

long double integrand(long double theta, long double a, long double b) {
    const long double s = std::sin(theta);
    const long double c = std::cos(theta);
    const long double aa = a * a;
    const long double bb = b * b;
    const long double e = std::sqrt((s * s) / aa + (c * c) / bb);

    const long double r = s * (a + 1.0L / (a * e));
    const long double minus_z_prime =
        s * (b + 1.0L / (b * e) + (c * c / b) * (1.0L / aa - 1.0L / bb) / (e * e * e));
    return r * r * minus_z_prime;
}

long double simpson(long double l, long double r, long double a, long double b) {
    const long double m = (l + r) * 0.5L;
    return (r - l) *
           (integrand(l, a, b) + 4.0L * integrand(m, a, b) + integrand(r, a, b)) / 6.0L;
}

long double adaptive_simpson(long double l,
                             long double r,
                             long double eps,
                             long double whole,
                             int depth,
                             long double a,
                             long double b) {
    const long double m = (l + r) * 0.5L;
    const long double left = simpson(l, m, a, b);
    const long double right = simpson(m, r, a, b);
    const long double delta = left + right - whole;

    if (depth <= 0 || std::fabsl(delta) <= 15.0L * eps) {
        return left + right + delta / 15.0L;
    }
    return adaptive_simpson(l, m, eps * 0.5L, left, depth - 1, a, b) +
           adaptive_simpson(m, r, eps * 0.5L, right, depth - 1, a, b);
}

long double compute_chocolate_volume(long double a, long double b) {
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const long double half_pi = pi * 0.5L;

    const long double whole = simpson(0.0L, half_pi, a, b);
    const long double integral = adaptive_simpson(0.0L, half_pi, 1e-14L, whole, 24, a, b);

    const long double outer = 2.0L * pi * integral;
    const long double inner = 4.0L * pi * a * a * b / 3.0L;
    return outer - inner;
}

bool run_checkpoints() {
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const long double c11 = compute_chocolate_volume(1.0L, 1.0L);
    const long double expected11 = 28.0L * pi / 3.0L;
    if (std::fabsl(c11 - expected11) > 1e-12L) {
        std::cerr << "Checkpoint failed: (a,b)=(1,1)\n";
        return false;
    }

    const long double c21 = compute_chocolate_volume(2.0L, 1.0L);
    if (std::fabsl(c21 - 60.35475635L) > 5e-9L) {
        std::cerr << "Checkpoint failed: (a,b)=(2,1)\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const long double answer = compute_chocolate_volume(options.a, options.b);
    std::cout << std::fixed << std::setprecision(8) << answer << '\n';
    return 0;
}
