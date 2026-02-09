#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

struct Vec2 {
    long double x{0.0L};
    long double y{0.0L};
};

Vec2 operator+(const Vec2& a, const Vec2& b) {
    return Vec2{a.x + b.x, a.y + b.y};
}

Vec2 operator-(const Vec2& a, const Vec2& b) {
    return Vec2{a.x - b.x, a.y - b.y};
}

Vec2 operator*(const Vec2& a, const long double k) {
    return Vec2{a.x * k, a.y * k};
}

long double dot(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}

long double norm(const Vec2& a) {
    return std::sqrt(a.x * a.x + a.y * a.y);
}

Vec2 i_times(const Vec2& z) {
    return Vec2{-z.y, z.x};
}

Vec2 complex_mul(const Vec2& z, const Vec2& c) {
    return Vec2{z.x * c.x - z.y * c.y, z.x * c.y + z.y * c.x};
}

// For c = u + iv, A_c = [[u,-v],[v,u]], so A_c^T d = [u dx + v dy, -v dx + u dy].
Vec2 apply_transpose_mul(const Vec2& d, const Vec2& c) {
    return Vec2{c.x * d.x + c.y * d.y, -c.y * d.x + c.x * d.y};
}

long double support_unit_square(const Vec2& d) {
    return std::max({0.0L, d.x, d.y, d.x + d.y});
}

struct Node {
    long double upper_bound;
    long double offset;
    Vec2 dir;
};

struct NodeCmp {
    bool operator()(const Node& a, const Node& b) const {
        return a.upper_bound < b.upper_bound;
    }
};

class SupportSolver {
   public:
    SupportSolver()
        : a_{16.0L / 25.0L, 12.0L / 25.0L},
          b_{9.0L / 25.0L, -12.0L / 25.0L},
          t1_{0.0L, 1.0L},
          t2_{16.0L / 25.0L, 37.0L / 25.0L} {}

    long double support(const Vec2& direction, const long double tol = 1e-15L) const {
        constexpr long double kRadialBound = 5.0L;  // Global |z| upper bound for the full tree.

        long double best = support_unit_square(direction);
        std::priority_queue<Node, std::vector<Node>, NodeCmp> pq;
        pq.push(Node{norm(direction) * kRadialBound, 0.0L, direction});

        while (!pq.empty() && pq.top().upper_bound > best + tol) {
            const Node cur = pq.top();
            pq.pop();

            best = std::max(best, cur.offset + support_unit_square(cur.dir));

            const Vec2 dir1 = apply_transpose_mul(cur.dir, a_);
            const long double off1 = cur.offset + dot(cur.dir, t1_);
            const long double ub1 = off1 + norm(dir1) * kRadialBound;
            if (ub1 > best + tol) {
                pq.push(Node{ub1, off1, dir1});
            }

            const Vec2 dir2 = apply_transpose_mul(cur.dir, b_);
            const long double off2 = cur.offset + dot(cur.dir, t2_);
            const long double ub2 = off2 + norm(dir2) * kRadialBound;
            if (ub2 > best + tol) {
                pq.push(Node{ub2, off2, dir2});
            }
        }

        return best;
    }

    Vec2 a() const { return a_; }
    Vec2 b() const { return b_; }
    Vec2 t1() const { return t1_; }
    Vec2 t2() const { return t2_; }

   private:
    Vec2 a_;
    Vec2 b_;
    Vec2 t1_;
    Vec2 t2_;
};

struct AffineSquare {
    Vec2 p;  // base start
    Vec2 z;  // base vector
};

std::tuple<long double, long double, long double, long double> finite_depth_bbox(const int depth) {
    const Vec2 a{16.0L / 25.0L, 12.0L / 25.0L};
    const Vec2 b{9.0L / 25.0L, -12.0L / 25.0L};
    std::vector<AffineSquare> current;
    current.push_back(AffineSquare{Vec2{0.0L, 0.0L}, Vec2{1.0L, 0.0L}});

    long double xmin = std::numeric_limits<long double>::infinity();
    long double xmax = -std::numeric_limits<long double>::infinity();
    long double ymin = std::numeric_limits<long double>::infinity();
    long double ymax = -std::numeric_limits<long double>::infinity();

    for (int level = 0; level <= depth; ++level) {
        for (const AffineSquare& sq : current) {
            const Vec2 iz = i_times(sq.z);
            const Vec2 c0 = sq.p;
            const Vec2 c1 = sq.p + sq.z;
            const Vec2 c2 = sq.p + iz;
            const Vec2 c3 = sq.p + sq.z + iz;

            for (const Vec2& c : {c0, c1, c2, c3}) {
                xmin = std::min(xmin, c.x);
                xmax = std::max(xmax, c.x);
                ymin = std::min(ymin, c.y);
                ymax = std::max(ymax, c.y);
            }
        }

        if (level == depth) {
            break;
        }

        std::vector<AffineSquare> next;
        next.reserve(current.size() * 2);
        for (const AffineSquare& sq : current) {
            const Vec2 iz = i_times(sq.z);
            const Vec2 left_p = sq.p + iz;
            const Vec2 left_z = complex_mul(sq.z, a);

            const Vec2 right_p = sq.p + iz + left_z;
            const Vec2 right_z = complex_mul(sq.z, b);

            next.push_back(AffineSquare{left_p, left_z});
            next.push_back(AffineSquare{right_p, right_z});
        }
        current.swap(next);
    }

    return {xmin, xmax, ymin, ymax};
}

bool nearly_equal(const long double a, const long double b, const long double tol = 1e-12L) {
    const long double scale = std::max(std::fabsl(a), std::fabsl(b));
    if (scale == 0.0L) {
        return true;
    }
    return std::fabsl(a - b) <= tol * scale;
}

bool run_checkpoints() {
    const auto [xmin1, xmax1, ymin1, ymax1] = finite_depth_bbox(1);
    if (!nearly_equal(xmin1, -12.0L / 25.0L) ||
        !nearly_equal(xmax1, 37.0L / 25.0L) ||
        !nearly_equal(ymin1, 0.0L) ||
        !nearly_equal(ymax1, 53.0L / 25.0L)) {
        std::cerr << "Checkpoint failed: finite depth-1 geometry\n";
        return false;
    }

    const SupportSolver solver;
    const std::vector<Vec2> dirs{{1.0L, 0.0L}, {-1.0L, 0.0L}, {0.0L, 1.0L}, {0.0L, -1.0L}};
    for (const Vec2& d : dirs) {
        const long double h = solver.support(d);
        const long double rhs = std::max({
            support_unit_square(d),
            dot(d, solver.t1()) + solver.support(apply_transpose_mul(d, solver.a())),
            dot(d, solver.t2()) + solver.support(apply_transpose_mul(d, solver.b())),
        });
        if (!nearly_equal(h, rhs, 2e-13L)) {
            std::cerr << "Checkpoint failed: support fixed-point residual\n";
            return false;
        }
    }

    return true;
}

long double solve_area() {
    const SupportSolver solver;
    const long double xmax = solver.support(Vec2{1.0L, 0.0L});
    const long double xmin = -solver.support(Vec2{-1.0L, 0.0L});
    const long double ymax = solver.support(Vec2{0.0L, 1.0L});
    const long double ymin = -solver.support(Vec2{0.0L, -1.0L});
    return (xmax - xmin) * (ymax - ymin);
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

    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(solve_area()) << '\n';
    return 0;
}
