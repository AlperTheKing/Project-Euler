#include <cassert>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>

namespace {

constexpr double kEps = 1e-12;

struct Point {
    double x;
    double y;
};

double distance(const Point a, const Point b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

Point circumcenter(const Point p1, const Point p2, const Point p3) {
    const double x1 = p1.x;
    const double y1 = p1.y;
    const double x2 = p2.x;
    const double y2 = p2.y;
    const double x3 = p3.x;
    const double y3 = p3.y;

    const double den = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

    const double n1 = x1 * x1 + y1 * y1;
    const double n2 = x2 * x2 + y2 * y2;
    const double n3 = x3 * x3 + y3 * y3;

    const double ux = (n1 * (y2 - y3) + n2 * (y3 - y1) + n3 * (y1 - y2)) / den;
    const double uy = (n1 * (x3 - x2) + n2 * (x1 - x3) + n3 * (x2 - x1)) / den;

    return {ux, uy};
}

double d_value(const int ra, const int rb, const int rc) {
    const double a = static_cast<double>(ra);
    const double b = static_cast<double>(rb);
    const double c = static_cast<double>(rc);

    const Point A{0.0, 0.0};
    const Point B{a + b, 0.0};

    const double cx = ((a + c) * (a + c) + (a + b) * (a + b) - (b + c) * (b + c)) / (2.0 * (a + b));
    const double cy = std::sqrt((a + c) * (a + c) - cx * cx);
    const Point C{cx, cy};

    const Point Pab{a, 0.0};
    const Point Pac{a / (a + c) * cx, a / (a + c) * cy};
    const Point Pbc{B.x + b / (b + c) * (cx - B.x), b / (b + c) * cy};

    const Point D = circumcenter(Pab, Pac, Pbc);

    const double k1 = 1.0 / a;
    const double k2 = 1.0 / b;
    const double k3 = 1.0 / c;
    const double k4 = k1 + k2 + k3 + 2.0 * std::sqrt(k1 * k2 + k2 * k3 + k3 * k1);
    const double r = 1.0 / k4;

    const double rhs_ab = (B.x * B.x + B.y * B.y) + (a + r) * (a + r) - (b + r) * (b + r);
    const double ex = rhs_ab / (2.0 * B.x);

    const double rhs_ac = (C.x * C.x + C.y * C.y) + (a + r) * (a + r) - (c + r) * (c + r);
    const double ey = (rhs_ac - 2.0 * C.x * ex) / (2.0 * C.y);
    const Point E{ex, ey};

    assert(std::fabs(distance(E, A) - (a + r)) < 1e-9);
    assert(std::fabs(distance(E, B) - (b + r)) < 1e-9);
    assert(std::fabs(distance(E, C) - (c + r)) < 1e-9);

    return distance(D, E);
}

double expected_d() {
    double sum = 0.0;
    std::int64_t count = 0;

    for (int a = 1; a <= 100; ++a) {
        for (int b = a + 1; b <= 100; ++b) {
            for (int c = b + 1; c <= 100; ++c) {
                if (std::gcd(a, std::gcd(b, c)) != 1) {
                    continue;
                }
                sum += d_value(a, b, c);
                ++count;
            }
        }
    }

    assert(count == 135739);
    return sum / static_cast<double>(count);
}

}  // namespace

int main() {
    assert(std::fabs(d_value(1, 2, 3) - 0.15676309893321677) < 1e-12);
    assert(std::fabs(d_value(1, 2, 4) - 0.19528865932762782) < 1e-12);
    assert(std::fabs(d_value(2, 3, 5) - 0.2060868349032223) < 1e-12);

    std::cout << std::fixed << std::setprecision(8) << expected_d() << '\n';
    return 0;
}
