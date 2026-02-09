#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int n = 400;
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
    int parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 2 && (options.n % 2 == 0);
}

long double f(const long double x) {
    const long double inside = 1.0L - x * x;
    return std::sqrt(std::max(0.0L, inside));
}

bool build_sequence(const int m, const long double x1, std::vector<long double>& x) {
    x.assign(static_cast<std::size_t>(m + 2), 0.0L);
    x[0] = 0.0L;
    x[1] = x1;

    if (!(x1 > 0.0L && x1 < 1.0L)) {
        return false;
    }

    for (int k = 1; k <= m; ++k) {
        const long double x_prev = x[static_cast<std::size_t>(k - 1)];
        const long double x_cur = x[static_cast<std::size_t>(k)];
        const long double y_prev = f(x_prev);
        const long double y_cur = f(x_cur);

        if (x_cur <= 0.0L || y_cur <= 0.0L) {
            return false;
        }

        const long double x_next = x_cur + (y_prev - y_cur) * y_cur / x_cur;
        if (!std::isfinite(static_cast<double>(x_next))) {
            return false;
        }
        x[static_cast<std::size_t>(k + 1)] = x_next;
    }
    return true;
}

long double shooting_residual(const int m, const long double x1) {
    std::vector<long double> x;
    if (!build_sequence(m, x1, x)) {
        return std::numeric_limits<long double>::infinity();
    }
    return x[static_cast<std::size_t>(m + 1)] - 1.0L;
}

std::vector<long double> solve_breakpoints(const int n) {
    const int m = n / 2;
    long double lo = 1e-18L;
    long double hi = 0.999999999999L;

    for (int iter = 0; iter < 260; ++iter) {
        const long double mid = (lo + hi) / 2.0L;
        const long double r = shooting_residual(m, mid);
        if (r > 0.0L) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    const long double x1 = (lo + hi) / 2.0L;
    std::vector<long double> x;
    build_sequence(m, x1, x);
    x[static_cast<std::size_t>(m + 1)] = 1.0L;  // enforce boundary exactly after shooting.
    return x;
}

long double minimal_red_area(const int n) {
    const int m = n / 2;
    const std::vector<long double> x = solve_breakpoints(n);

    long double quadrant_area = 0.0L;
    for (int i = 0; i <= m; ++i) {
        const long double width = x[static_cast<std::size_t>(i + 1)] - x[static_cast<std::size_t>(i)];
        quadrant_area += width * f(x[static_cast<std::size_t>(i)]);
    }

    return 4.0L * quadrant_area;
}

bool nearly_equal(long double a, long double b, long double tol = 1e-13L) {
    const long double scale = std::max(std::fabsl(a), std::fabsl(b));
    if (scale == 0.0L) {
        return true;
    }
    return std::fabsl(a - b) <= tol * scale;
}

bool run_checkpoints() {
    const long double area_n2 = minimal_red_area(2);
    const long double expected_n2 = 4.0L * std::sqrt(2.0L) - 2.0L;
    if (!nearly_equal(area_n2, expected_n2, 1e-12L)) {
        std::cerr << "Checkpoint failed for N=2\n";
        return false;
    }

    const long double area_n10 = minimal_red_area(10);
    const long double expected_n10 = 3.3469640797L;
    if (!nearly_equal(area_n10, expected_n10, 2e-11L)) {
        std::cerr << "Checkpoint failed for N=10 (got " << std::setprecision(16)
                  << static_cast<double>(area_n10) << ")\n";
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

    const long double answer = minimal_red_area(options.n);
    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(answer) << '\n';
    return 0;
}
