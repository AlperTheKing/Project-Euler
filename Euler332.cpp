#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;

struct Point {
    int x;
    int y;
    int z;
};

std::vector<Point> lattice_points_on_sphere(const int r) {
    const int rr = r * r;
    std::vector<Point> points;
    points.reserve(600);

    for (int x = -r; x <= r; ++x) {
        const int x2 = x * x;
        for (int y = -r; y <= r; ++y) {
            const int y2 = y * y;
            const int z2 = rr - x2 - y2;
            if (z2 < 0) {
                continue;
            }
            const int z = static_cast<int>(std::llround(std::sqrt(static_cast<long double>(z2))));
            if (z * z != z2) {
                continue;
            }
            points.push_back({x, y, z});
            if (z != 0) {
                points.push_back({x, y, -z});
            }
        }
    }

    return points;
}

i64 dot(const Point& a, const Point& b) {
    return static_cast<i64>(a.x) * b.x + static_cast<i64>(a.y) * b.y + static_cast<i64>(a.z) * b.z;
}

i64 abs_det(const Point& a, const Point& b, const Point& c) {
    const i64 det = static_cast<i64>(a.x) * (static_cast<i64>(b.y) * c.z - static_cast<i64>(b.z) * c.y) -
                    static_cast<i64>(a.y) * (static_cast<i64>(b.x) * c.z - static_cast<i64>(b.z) * c.x) +
                    static_cast<i64>(a.z) * (static_cast<i64>(b.x) * c.y - static_cast<i64>(b.y) * c.x);
    return det >= 0 ? det : -det;
}

long double min_spherical_triangle_area(const int r) {
    const std::vector<Point> points = lattice_points_on_sphere(r);
    const int n = static_cast<int>(points.size());
    if (n < 3) {
        return 0.0L;
    }

    std::vector<std::vector<i64>> dots(static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const i64 v = dot(points[static_cast<std::size_t>(i)], points[static_cast<std::size_t>(j)]);
            dots[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = v;
            dots[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = v;
        }
    }

    i64 best_det = std::numeric_limits<i64>::max();
    i64 best_den = 1;

    const i64 r2 = static_cast<i64>(r) * r;
    const i64 r3 = r2 * static_cast<i64>(r);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                const i64 detv = abs_det(points[static_cast<std::size_t>(i)], points[static_cast<std::size_t>(j)],
                                         points[static_cast<std::size_t>(k)]);
                if (detv == 0) {
                    continue;
                }

                const i64 sum_dots = dots[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                                     dots[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] +
                                     dots[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
                const i64 denv = r3 + static_cast<i64>(r) * sum_dots;

                if (denv <= 0) {
                    continue;
                }

                // Minimize det/den for positive numerator/denominator.
                if (best_det == std::numeric_limits<i64>::max() ||
                    static_cast<i128>(detv) * static_cast<i128>(best_den) <
                        static_cast<i128>(best_det) * static_cast<i128>(denv)) {
                    best_det = detv;
                    best_den = denv;
                }
            }
        }
    }

    if (best_det == std::numeric_limits<i64>::max()) {
        return 0.0L;
    }

    const long double omega = 2.0L * std::atan2(static_cast<long double>(best_det), static_cast<long double>(best_den));
    return static_cast<long double>(r2) * omega;
}

long double solve_sum() {
    long double total = 0.0L;
    for (int r = 1; r <= 50; ++r) {
        total += min_spherical_triangle_area(r);
    }
    return total;
}

bool run_checkpoints() {
    const long double a14 = min_spherical_triangle_area(14);
    const long double expected = 3.294040L;
    if (std::fabsl(a14 - expected) > 0.0000005L) {
        std::cerr << "Checkpoint failed: A(14) mismatch, got=" << std::fixed << std::setprecision(9) << a14 << '\n';
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

    const long double answer = solve_sum();
    std::cout << std::fixed << std::setprecision(6) << answer << '\n';
    return 0;
}
