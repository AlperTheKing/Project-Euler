#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

long double p4(long double x) {
    const long double x2 = x * x;
    return x2 * x2;
}

std::vector<long double> optimal_knots(int n) {
    const int m = n - 1;
    std::vector<long double> x(m + 1, 0.0L);
    for (int i = 0; i <= m; ++i) {
        x[i] = -1.0L + 2.0L * static_cast<long double>(i) / static_cast<long double>(m);
    }

    constexpr long double omega = 0.5L;
    constexpr long double tol = 1e-18L;
    constexpr int max_iter = 300000;

    for (int it = 0; it < max_iter; ++it) {
        long double max_delta = 0.0L;
        for (int i = 1; i < m; ++i) {
            const long double a = x[i - 1];
            const long double b = x[i + 1];
            const long double target = std::cbrt((p4(b) - p4(a)) / (4.0L * (b - a)));
            const long double nx = x[i] + omega * (target - x[i]);
            const long double d = std::fabsl(nx - x[i]);
            if (d > max_delta) {
                max_delta = d;
            }
            x[i] = nx;
        }
        if (max_delta < tol) {
            break;
        }
    }

    for (int i = 0; i < m; ++i) {
        assert(x[i] < x[i + 1]);
    }
    long double max_residual = 0.0L;
    for (int i = 1; i < m; ++i) {
        const long double lhs = 4.0L * x[i] * x[i] * x[i] * (x[i + 1] - x[i - 1]);
        const long double rhs = p4(x[i + 1]) - p4(x[i - 1]);
        max_residual = std::max(max_residual, std::fabsl(lhs - rhs));
    }
    assert(max_residual < 1e-15L);

    return x;
}

long double polygon_area(const std::vector<long double>& x) {
    const int m = static_cast<int>(x.size()) - 1;
    long double area = 0.0L;
    for (int i = 0; i < m; ++i) {
        const long double dx = x[i + 1] - x[i];
        area += (2.0L - (p4(x[i]) + p4(x[i + 1]))) * dx * 0.5L;
    }
    return area;
}

long double G(int n) {
    return polygon_area(optimal_knots(n));
}

void validate() {
    const long double g3 = G(3);
    assert(std::fabsl(g3 - 1.0L) < 1e-15L);

    const long double g5 = G(5);
    assert(std::fabsl(g5 - 1.477309771L) < 1e-9L);
}

}  // namespace

int main() {
    validate();
    std::cout << std::fixed << std::setprecision(9) << G(101) << '\n';
    return 0;
}
