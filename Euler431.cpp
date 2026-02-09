#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

constexpr long double PI = 3.141592653589793238462643383279502884L;

struct Options {
    long double radius = 6.0L;
    long double alpha_deg = 40.0L;
    int simpson_n = 4096;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_ld_after_prefix(const std::string& arg, const std::string& prefix, long double& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stold(tail);
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
        if (parse_ld_after_prefix(arg, "--radius=", options.radius) ||
            parse_ld_after_prefix(arg, "--alpha=", options.alpha_deg) ||
            parse_int_after_prefix(arg, "--simpson-n=", options.simpson_n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.radius > 0.0L && options.alpha_deg > 0.0L && options.alpha_deg < 90.0L &&
           options.simpson_n >= 128 && (options.simpson_n % 2 == 0);
}

long double integral_distance(long double x, long double radius, int simpson_n) {
    // I(x) = ∬_disk dist((u,v),(x,0)) dudv
    //      = (1/3) ∫_0^{2π} s(phi)^3 dphi,
    // where s(phi) is radial intersection in polar coordinates centered at (x,0).
    const long double a = 0.0L;
    const long double b = 2.0L * PI;
    const long double h = (b - a) / static_cast<long double>(simpson_n);

    long double acc = 0.0L;
    for (int i = 0; i <= simpson_n; ++i) {
        const long double phi = a + h * static_cast<long double>(i);
        const long double sinp = std::sin(phi);
        const long double cosp = std::cos(phi);
        const long double rad = std::sqrt(std::max(0.0L, radius * radius - x * x * sinp * sinp));
        const long double s = -x * cosp + rad;
        const long double f = s * s * s;

        int coeff = 2;
        if (i == 0 || i == simpson_n) {
            coeff = 1;
        } else if (i & 1) {
            coeff = 4;
        }
        acc += static_cast<long double>(coeff) * f;
    }

    return (h / 3.0L) * acc / 3.0L;
}

long double wasted_volume(long double x, long double radius, long double alpha_deg, int simpson_n) {
    const long double tan_alpha = std::tan(alpha_deg * PI / 180.0L);
    return tan_alpha * integral_distance(x, radius, simpson_n);
}

std::vector<long double> solve_x_values(long double radius, long double alpha_deg, int simpson_n) {
    const long double v0 = wasted_volume(0.0L, radius, alpha_deg, simpson_n);
    const long double vr = wasted_volume(radius, radius, alpha_deg, simpson_n);

    const int k_lo = static_cast<int>(std::ceill(std::sqrt(v0) - 1e-14L));
    const int k_hi = static_cast<int>(std::floorl(std::sqrt(vr) + 1e-14L));

    std::vector<long double> xs;
    for (int k = k_lo; k <= k_hi; ++k) {
        const long double target = static_cast<long double>(k) * static_cast<long double>(k);
        if (target < v0 - 1e-12L || target > vr + 1e-12L) {
            continue;
        }

        long double lo = 0.0L;
        long double hi = radius;
        for (int it = 0; it < 120; ++it) {
            const long double mid = (lo + hi) * 0.5L;
            const long double vm = wasted_volume(mid, radius, alpha_deg, simpson_n);
            if (vm < target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        xs.push_back((lo + hi) * 0.5L);
    }

    return xs;
}

bool close_to(long double a, long double b, long double tol) {
    return std::fabsl(a - b) <= tol;
}

bool run_checkpoints(const int simpson_n) {
    const long double v0 = wasted_volume(0.0L, 3.0L, 30.0L, simpson_n);
    if (!close_to(v0, 32.648388556L, 5e-10L)) {
        std::cerr << "Checkpoint failed: V(0) for R=3, alpha=30\n";
        return false;
    }

    const auto xs = solve_x_values(3.0L, 30.0L, simpson_n);
    if (xs.size() != 2U) {
        std::cerr << "Checkpoint failed: expected two square solutions for sample case\n";
        return false;
    }
    if (!close_to(xs[0], 1.114785284L, 5e-10L)) {
        std::cerr << "Checkpoint failed: sample x for V=36\n";
        return false;
    }
    if (!close_to(xs[1], 2.511167869L, 5e-10L)) {
        std::cerr << "Checkpoint failed: sample x for V=49\n";
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
    if (options.run_checkpoints && !run_checkpoints(options.simpson_n)) {
        return 2;
    }

    const auto xs = solve_x_values(options.radius, options.alpha_deg, options.simpson_n);
    long double sum = 0.0L;
    for (long double x : xs) {
        sum += x;
    }

    std::cout << std::fixed << std::setprecision(9) << static_cast<double>(sum) << '\n';
    return 0;
}
