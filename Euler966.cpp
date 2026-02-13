#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

namespace {

constexpr long double kPi = 3.141592653589793238462643383279502884L;
constexpr long double kEps = 1e-12L;

struct Vec2 {
    long double x;
    long double y;

    Vec2 operator+(const Vec2& other) const {
        return {x + other.x, y + other.y};
    }

    Vec2 operator-(const Vec2& other) const {
        return {x - other.x, y - other.y};
    }

    Vec2 operator*(long double s) const {
        return {x * s, y * s};
    }
};

long double dot(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}

long double cross(const Vec2& a, const Vec2& b) {
    return a.x * b.y - a.y * b.x;
}

long double norm2(const Vec2& v) {
    return dot(v, v);
}

long double dist(const Vec2& a, const Vec2& b) {
    return std::sqrt(norm2(a - b));
}

struct Triangle {
    Vec2 a;
    Vec2 b;
    Vec2 c;
    long double area;
    long double radius;
    long double max_side;
};

Triangle make_triangle(int s1, int s2, int s3) {
    std::array<int, 3> sides{s1, s2, s3};
    std::sort(sides.begin(), sides.end());

    const long double a = static_cast<long double>(sides[0]);
    const long double b = static_cast<long double>(sides[1]);
    const long double c = static_cast<long double>(sides[2]);

    const Vec2 A{0.0L, 0.0L};
    const Vec2 B{c, 0.0L};
    const long double x = (b * b + c * c - a * a) / (2.0L * c);
    long double y2 = b * b - x * x;
    if (y2 < 0.0L && y2 > -1e-18L) {
        y2 = 0.0L;
    }
    const Vec2 C{x, std::sqrt(std::max(0.0L, y2))};

    const long double area = std::abs(cross(B - A, C - A)) * 0.5L;
    const long double radius = std::sqrt(area / kPi);
    const long double max_side = std::max({dist(A, B), dist(B, C), dist(C, A)});
    return {A, B, C, area, radius, max_side};
}

bool point_in_triangle(const Triangle& tri, const Vec2& p) {
    const long double c1 = cross(tri.b - tri.a, p - tri.a);
    const long double c2 = cross(tri.c - tri.b, p - tri.b);
    const long double c3 = cross(tri.a - tri.c, p - tri.c);
    return c1 >= -1e-14L && c2 >= -1e-14L && c3 >= -1e-14L;
}

Vec2 nearest_on_segment(const Vec2& p, const Vec2& u, const Vec2& v) {
    const Vec2 uv = v - u;
    const long double d = norm2(uv);
    if (d <= 0.0L) {
        return u;
    }
    long double t = dot(p - u, uv) / d;
    if (t < 0.0L) {
        t = 0.0L;
    } else if (t > 1.0L) {
        t = 1.0L;
    }
    return u + uv * t;
}

Vec2 project_to_triangle(const Triangle& tri, const Vec2& p) {
    if (point_in_triangle(tri, p)) {
        return p;
    }

    const Vec2 candidates[6] = {
        tri.a,
        tri.b,
        tri.c,
        nearest_on_segment(p, tri.a, tri.b),
        nearest_on_segment(p, tri.b, tri.c),
        nearest_on_segment(p, tri.c, tri.a),
    };

    Vec2 best = candidates[0];
    long double best_d = norm2(p - best);
    for (int i = 1; i < 6; ++i) {
        const long double cur = norm2(p - candidates[i]);
        if (cur < best_d) {
            best_d = cur;
            best = candidates[i];
        }
    }
    return best;
}

long double segment_contribution(const Vec2& a, const Vec2& b, long double r2) {
    Vec2 pts[4];
    int n = 0;
    pts[n++] = a;

    const Vec2 d = b - a;
    const long double qa = dot(d, d);
    const long double qb = 2.0L * dot(a, d);
    const long double qc = dot(a, a) - r2;
    const long double disc = qb * qb - 4.0L * qa * qc;

    if (disc > 1e-18L) {
        const long double sdisc = std::sqrt(disc);
        long double t1 = (-qb - sdisc) / (2.0L * qa);
        long double t2 = (-qb + sdisc) / (2.0L * qa);
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        if (t1 > kEps && t1 < 1.0L - kEps) {
            pts[n++] = a + d * t1;
        }
        if (t2 > kEps && t2 < 1.0L - kEps && std::abs(t2 - t1) > 1e-12L) {
            pts[n++] = a + d * t2;
        }
    }

    pts[n++] = b;

    long double acc = 0.0L;
    for (int i = 0; i + 1 < n; ++i) {
        const Vec2 p = pts[i];
        const Vec2 q = pts[i + 1];
        const Vec2 mid{(p.x + q.x) * 0.5L, (p.y + q.y) * 0.5L};
        if (norm2(mid) <= r2 + 1e-13L) {
            acc += 0.5L * cross(p, q);
        } else {
            acc += 0.5L * r2 * std::atan2(cross(p, q), dot(p, q));
        }
    }
    return acc;
}

long double overlap_area(const Triangle& tri, const Vec2& center) {
    const long double r = tri.radius;
    const long double r2 = r * r;

    const Vec2 p0 = tri.a - center;
    const Vec2 p1 = tri.b - center;
    const Vec2 p2 = tri.c - center;

    long double area = 0.0L;
    area += segment_contribution(p0, p1, r2);
    area += segment_contribution(p1, p2, r2);
    area += segment_contribution(p2, p0, r2);
    return std::abs(area);
}

struct EvalPoint {
    Vec2 p;
    long double v;
};

EvalPoint nelder_mead(const Triangle& tri, const Vec2& start, long double step, int iterations) {
    auto eval = [&](const Vec2& raw) -> EvalPoint {
        const Vec2 projected = project_to_triangle(tri, raw);
        return {projected, overlap_area(tri, projected)};
    };

    std::array<EvalPoint, 3> simplex = {
        eval(start),
        eval({start.x + step, start.y}),
        eval({start.x, start.y + step}),
    };

    for (int it = 0; it < iterations; ++it) {
        std::sort(simplex.begin(), simplex.end(),
                  [](const EvalPoint& lhs, const EvalPoint& rhs) { return lhs.v > rhs.v; });

        const Vec2 centroid{
            (simplex[0].p.x + simplex[1].p.x) * 0.5L,
            (simplex[0].p.y + simplex[1].p.y) * 0.5L,
        };

        const Vec2 reflected{
            centroid.x + (centroid.x - simplex[2].p.x),
            centroid.y + (centroid.y - simplex[2].p.y),
        };
        const EvalPoint er = eval(reflected);

        if (er.v > simplex[0].v) {
            const Vec2 expanded{
                centroid.x + 2.0L * (er.p.x - centroid.x),
                centroid.y + 2.0L * (er.p.y - centroid.y),
            };
            const EvalPoint ee = eval(expanded);
            simplex[2] = (ee.v > er.v ? ee : er);
            continue;
        }

        if (er.v > simplex[1].v) {
            simplex[2] = er;
            continue;
        }

        const Vec2 contracted{
            centroid.x + 0.5L * (simplex[2].p.x - centroid.x),
            centroid.y + 0.5L * (simplex[2].p.y - centroid.y),
        };
        const EvalPoint ec = eval(contracted);
        if (ec.v > simplex[2].v) {
            simplex[2] = ec;
            continue;
        }

        simplex[1] = eval({
            simplex[0].p.x + 0.5L * (simplex[1].p.x - simplex[0].p.x),
            simplex[0].p.y + 0.5L * (simplex[1].p.y - simplex[0].p.y),
        });
        simplex[2] = eval({
            simplex[0].p.x + 0.5L * (simplex[2].p.x - simplex[0].p.x),
            simplex[0].p.y + 0.5L * (simplex[2].p.y - simplex[0].p.y),
        });
    }

    std::sort(simplex.begin(), simplex.end(),
              [](const EvalPoint& lhs, const EvalPoint& rhs) { return lhs.v > rhs.v; });
    return simplex[0];
}

long double maximize_overlap(int s1, int s2, int s3) {
    const Triangle tri = make_triangle(s1, s2, s3);

    std::vector<Vec2> seeds;
    seeds.reserve(24);

    const int grid = 5;
    for (int i = 0; i <= grid; ++i) {
        for (int j = 0; j + i <= grid; ++j) {
            const long double u = static_cast<long double>(i) / static_cast<long double>(grid);
            const long double v = static_cast<long double>(j) / static_cast<long double>(grid);
            const long double w = 1.0L - u - v;
            seeds.push_back({
                u * tri.a.x + v * tri.b.x + w * tri.c.x,
                u * tri.a.y + v * tri.b.y + w * tri.c.y,
            });
        }
    }

    const Vec2 centroid{
        (tri.a.x + tri.b.x + tri.c.x) / 3.0L,
        (tri.a.y + tri.b.y + tri.c.y) / 3.0L,
    };
    seeds.push_back(centroid);

    const long double side_a = dist(tri.b, tri.c);
    const long double side_b = dist(tri.c, tri.a);
    const long double side_c = dist(tri.a, tri.b);
    const long double ws = side_a + side_b + side_c;
    seeds.push_back({
        (side_a * tri.a.x + side_b * tri.b.x + side_c * tri.c.x) / ws,
        (side_a * tri.a.y + side_b * tri.b.y + side_c * tri.c.y) / ws,
    });

    std::vector<EvalPoint> evaluated;
    evaluated.reserve(seeds.size());
    for (const Vec2& p : seeds) {
        const Vec2 projected = project_to_triangle(tri, p);
        evaluated.push_back({projected, overlap_area(tri, projected)});
    }

    std::sort(evaluated.begin(), evaluated.end(),
              [](const EvalPoint& lhs, const EvalPoint& rhs) { return lhs.v > rhs.v; });

    EvalPoint best = evaluated[0];
    const int tries = std::min(4, static_cast<int>(evaluated.size()));
    const long double step = tri.max_side * 0.18L;
    for (int i = 0; i < tries; ++i) {
        const EvalPoint cur = nelder_mead(tri, evaluated[static_cast<std::size_t>(i)].p, step, 90);
        if (cur.v > best.v) {
            best = cur;
        }
    }

    const EvalPoint refined = nelder_mead(tri, best.p, tri.max_side * 0.05L, 50);
    if (refined.v > best.v) {
        best = refined;
    }

    if (best.v < 0.0L) {
        return 0.0L;
    }
    if (best.v > tri.area) {
        return tri.area;
    }
    return best.v;
}

bool nearly_equal(long double x, long double y, long double eps) {
    return std::abs(x - y) <= eps;
}

bool run_checkpoints() {
    const long double i345 = maximize_overlap(3, 4, 5);
    if (!nearly_equal(i345, 4.593049L, 2e-6L)) {
        std::cerr << "Checkpoint failed for I(3,4,5): got " << std::setprecision(12) << i345 << '\n';
        return false;
    }

    const long double i346 = maximize_overlap(3, 4, 6);
    if (!nearly_equal(i346, 3.552564L, 2e-6L)) {
        std::cerr << "Checkpoint failed for I(3,4,6): got " << std::setprecision(12) << i346 << '\n';
        return false;
    }

    const long double i_perm = maximize_overlap(5, 3, 4);
    if (!nearly_equal(i_perm, i345, 1e-10L)) {
        std::cerr << "Permutation checkpoint failed for I(5,3,4)" << '\n';
        return false;
    }

    return true;
}

std::vector<std::array<int, 3>> enumerate_triangles() {
    std::vector<std::array<int, 3>> all;
    for (int a = 1; a <= 200; ++a) {
        for (int b = a; b <= 200 - a; ++b) {
            const int max_c = 200 - a - b;
            for (int c = b; c <= max_c; ++c) {
                if (a + b <= c) {
                    continue;
                }
                all.push_back({a, b, c});
            }
        }
    }
    return all;
}

long double solve() {
    const std::vector<std::array<int, 3>> triangles = enumerate_triangles();

    unsigned int thread_count = std::thread::hardware_concurrency();
    if (thread_count == 0) {
        thread_count = 4;
    }
    thread_count = std::min<unsigned int>(thread_count, static_cast<unsigned int>(triangles.size()));
    if (thread_count == 0) {
        return 0.0L;
    }

    std::atomic<std::size_t> next_idx{0};
    std::vector<long double> partial(thread_count, 0.0L);

    auto worker = [&](unsigned int tid) {
        long double sum = 0.0L;
        long double comp = 0.0L;

        while (true) {
            const std::size_t idx = next_idx.fetch_add(1, std::memory_order_relaxed);
            if (idx >= triangles.size()) {
                break;
            }
            const auto& t = triangles[idx];
            const long double value = maximize_overlap(t[0], t[1], t[2]);

            const long double y = value - comp;
            const long double z = sum + y;
            comp = (z - sum) - y;
            sum = z;
        }

        partial[tid] = sum;
    };

    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (unsigned int tid = 0; tid < thread_count; ++tid) {
        threads.emplace_back(worker, tid);
    }
    for (std::thread& th : threads) {
        th.join();
    }

    long double total = 0.0L;
    long double comp = 0.0L;
    for (long double part : partial) {
        const long double y = part - comp;
        const long double z = total + y;
        comp = (z - total) - y;
        total = z;
    }
    return total;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    const long double answer = solve();
    std::cout << std::fixed << std::setprecision(2) << std::round(answer * 100.0L) / 100.0L << '\n';
    return 0;
}
