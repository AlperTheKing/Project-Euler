#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

constexpr double kEps = 1e-10;

struct Options {
    int size = 36;
    bool run_checkpoints = true;
};

struct Point {
    double x = 0.0;
    double y = 0.0;
};

struct Line {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    Line() = default;

    Line(const Point& p, const Point& q)
        : a(p.y - q.y), b(q.x - p.x), c(q.x * p.y - p.x * q.y) {}
};

struct Intersection {
    bool valid = false;
    Point p{};
};

bool almost_equal(double x, double y) {
    return std::fabs(x - y) <= kEps;
}

bool same_point(const Point& p, const Point& q) {
    return almost_equal(p.x, q.x) && almost_equal(p.y, q.y);
}

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_int_after_prefix(arg, "--size=", options.size)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.size >= 1;
}

Intersection intersect(const Line& lhs, const Line& rhs) {
    const double det = lhs.a * rhs.b - rhs.a * lhs.b;
    if (std::fabs(det) < kEps) {
        return {false, {}};
    }

    const double x = (lhs.c * rhs.b - rhs.c * lhs.b) / det;
    const double y = (lhs.a * rhs.c - rhs.a * lhs.c) / det;
    return {true, {x, y}};
}

double line_det(const Point& p, const Point& q, const Point& r) {
    return (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
}

bool inside_hull(const Point& p, const Point& a, const Point& b, const Point& c) {
    if (line_det(a, b, p) < -kEps) {
        return false;
    }
    if (line_det(b, c, p) < -kEps) {
        return false;
    }
    if (line_det(c, a, p) < -kEps) {
        return false;
    }
    return true;
}

std::vector<Line> build_lines(const int size) {
    const double h = std::sqrt(3.0) / 2.0;

    const Point a{0.0, 0.0};
    const Point b{1.0, 0.0};
    const Point c{0.5, h};

    const Point ab{(a.x + b.x) * 0.5, (a.y + b.y) * 0.5};
    const Point ac{(a.x + c.x) * 0.5, (a.y + c.y) * 0.5};
    const Point bc{(b.x + c.x) * 0.5, (b.y + c.y) * 0.5};

    std::vector<Line> lines;
    lines.reserve(static_cast<std::size_t>(9 * size - 3));

    // Family parallel to AB.
    for (int i = 0; i < size; ++i) {
        const double y = static_cast<double>(i) * h;
        lines.emplace_back(Point{a.x, y}, Point{b.x, y});
    }

    // Family through A and BC-like directions (both sides).
    for (int i = 0; i < size; ++i) {
        const double x = static_cast<double>(i);
        lines.emplace_back(Point{x, a.y}, Point{bc.x + x, bc.y});
        if (i > 0) {
            lines.emplace_back(Point{-x, a.y}, Point{bc.x - x, bc.y});
        }
    }

    // Family parallel to AC.
    for (int i = 0; i < size; ++i) {
        const double x = static_cast<double>(i);
        lines.emplace_back(Point{x, a.y}, Point{c.x + x, c.y});
    }

    // Family parallel to BC.
    for (int i = 0; i < size; ++i) {
        const double x = static_cast<double>(i);
        lines.emplace_back(Point{x + 1.0, b.y}, Point{c.x + x, c.y});
    }

    // Family through B and AC-like directions.
    for (int i = 0; i < 2 * size - 1; ++i) {
        const double x = static_cast<double>(i);
        lines.emplace_back(Point{x + 1.0, b.y}, Point{ac.x + x, ac.y});
    }

    // Family through C and AB-like directions.
    for (int i = 1; i < 2 * size; ++i) {
        const double x = static_cast<double>(i) * c.x;
        lines.emplace_back(Point{x, 0.0}, Point{x, h});
    }

    return lines;
}

i64 solve(const int size) {
    const double h = std::sqrt(3.0) / 2.0;
    const Point a{0.0, 0.0};
    const Point b{static_cast<double>(size), 0.0};
    const Point c{0.5 * static_cast<double>(size), h * static_cast<double>(size)};

    const std::vector<Line> lines = build_lines(size);
    const int n = static_cast<int>(lines.size());

    std::vector<Intersection> pair_data(static_cast<std::size_t>(n * n));
    auto at = [&](int i, int j) -> Intersection& {
        return pair_data[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                         static_cast<std::size_t>(j)];
    };

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Intersection inter = intersect(lines[static_cast<std::size_t>(i)],
                                           lines[static_cast<std::size_t>(j)]);
            if (inter.valid && !inside_hull(inter.p, a, b, c)) {
                inter.valid = false;
            }
            at(i, j) = inter;
        }
    }

    i64 count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const Intersection& ij = at(i, j);
            if (!ij.valid) {
                continue;
            }
            for (int k = j + 1; k < n; ++k) {
                const Intersection& ik = at(i, k);
                const Intersection& jk = at(j, k);
                if (!ik.valid || !jk.valid) {
                    continue;
                }
                if (same_point(ij.p, ik.p) || same_point(ij.p, jk.p) || same_point(ik.p, jk.p)) {
                    continue;
                }
                ++count;
            }
        }
    }

    return count;
}

bool run_checkpoints() {
    if (solve(1) != 16) {
        std::cerr << "Checkpoint failed for size=1" << '\n';
        return false;
    }
    if (solve(2) != 104) {
        std::cerr << "Checkpoint failed for size=2" << '\n';
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
    std::cout << solve(options.size) << '\n';
    return 0;
}
