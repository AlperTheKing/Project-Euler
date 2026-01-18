#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

namespace {

struct Point {
    long double x;
    long double y;
};

constexpr long double kPi = 3.141592653589793238462643383279502884L;
constexpr long double kEps = 1e-13L;
constexpr long double kExpected = 1.5276527928L;

long double dist(const Point& a, const Point& b) {
    long double dx = a.x - b.x;
    long double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

long double polygon_area(const std::vector<Point>& poly) {
    long double s = 0.0L;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const Point& p = poly[i];
        const Point& q = poly[(i + 1) % poly.size()];
        s += p.x * q.y - q.x * p.y;
    }
    return std::fabsl(s) * 0.5L;
}

long double orient(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool on_segment(const Point& a, const Point& b, const Point& c) {
    return std::min(a.x, b.x) - kEps <= c.x && c.x <= std::max(a.x, b.x) + kEps &&
           std::min(a.y, b.y) - kEps <= c.y && c.y <= std::max(a.y, b.y) + kEps;
}

bool segments_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    long double o1 = orient(a, b, c);
    long double o2 = orient(a, b, d);
    long double o3 = orient(c, d, a);
    long double o4 = orient(c, d, b);
    if (std::fabsl(o1) < kEps && on_segment(a, b, c)) return true;
    if (std::fabsl(o2) < kEps && on_segment(a, b, d)) return true;
    if (std::fabsl(o3) < kEps && on_segment(c, d, a)) return true;
    if (std::fabsl(o4) < kEps && on_segment(c, d, b)) return true;
    return (o1 * o2 < 0) && (o3 * o4 < 0);
}

bool is_simple_polygon(const std::vector<Point>& poly) {
    std::size_t n = poly.size();
    for (std::size_t i = 0; i < n; ++i) {
        Point a1 = poly[i];
        Point a2 = poly[(i + 1) % n];
        for (std::size_t j = i + 1; j < n; ++j) {
            if (j == i || (j + 1) % n == i || j == (i + 1) % n) {
                continue;
            }
            Point b1 = poly[j];
            Point b2 = poly[(j + 1) % n];
            if (segments_intersect(a1, a2, b1, b2)) {
                return false;
            }
        }
    }
    return true;
}

std::vector<std::vector<Point>> build_candidates(long double a) {
    long double b = kPi - a;
    Point v0{0.0L, 0.0L};
    Point v1{1.0L, 0.0L};
    Point v2{1.0L + std::cosl(a), std::sinl(a)};
    Point v4{std::cosl(b), std::sinl(b)};

    long double dx = v4.x - v2.x;
    long double dy = v4.y - v2.y;
    long double d = std::sqrt(dx * dx + dy * dy);
    long double half = 0.5L * d;
    long double h = std::sqrt(std::max(0.0L, 1.0L - half * half));
    long double mx = (v2.x + v4.x) * 0.5L;
    long double my = (v2.y + v4.y) * 0.5L;
    long double ox = -dy / d * h;
    long double oy = dx / d * h;

    Point p1{mx + ox, my + oy};
    Point p2{mx - ox, my - oy};
    std::vector<Point> poly1{v0, v1, v2, p1, v4};
    std::vector<Point> poly2{v0, v1, v2, p2, v4};
    return {poly1, poly2};
}

struct WidthData {
    std::vector<long double> w;
    std::vector<long double> wx;
    std::vector<long double> wy;
};

WidthData compute_widths(const std::vector<Point>& poly,
                         const std::vector<long double>& cos_t,
                         const std::vector<long double>& sin_t) {
    std::size_t n = cos_t.size();
    WidthData data;
    data.w.assign(n, 0.0L);
    data.wx.assign(n, 0.0L);
    data.wy.assign(n, 0.0L);

    std::size_t threads = std::max<unsigned>(1U, std::thread::hardware_concurrency());
    std::size_t chunk = (n + threads - 1) / threads;
    std::vector<std::thread> workers;
    workers.reserve(threads);

    auto worker = [&](std::size_t start) {
        std::size_t end = std::min(n, start + chunk);
        std::vector<Point> pts(poly.size());
        std::vector<Point> shifted(poly.size());
        for (std::size_t i = start; i < end; ++i) {
            long double c = cos_t[i];
            long double s = sin_t[i];
            long double minx = 1e100L, miny = 1e100L;
            long double maxx = -1e100L, maxy = -1e100L;
            for (std::size_t k = 0; k < poly.size(); ++k) {
                long double x = c * poly[k].x - s * poly[k].y;
                long double y = s * poly[k].x + c * poly[k].y;
                pts[k] = {x, y};
                minx = std::min(minx, x);
                miny = std::min(miny, y);
                maxx = std::max(maxx, x);
                maxy = std::max(maxy, y);
            }
            for (std::size_t k = 0; k < poly.size(); ++k) {
                shifted[k] = {pts[k].x - minx, pts[k].y - miny};
            }

            long double maxv = 0.0L;
            for (const auto& p : shifted) {
                long double v = (p.x < p.y) ? p.x : p.y;
                if (v > maxv) {
                    maxv = v;
                }
            }
            for (std::size_t k = 0; k < shifted.size(); ++k) {
                const Point& a = shifted[k];
                const Point& b = shifted[(k + 1) % shifted.size()];
                long double d1 = a.x - a.y;
                long double d2 = b.x - b.y;
                if (std::fabsl(d1) < kEps && std::fabsl(d2) < kEps) {
                    maxv = std::max(maxv, a.x);
                } else if (std::fabsl(d1) < kEps) {
                    maxv = std::max(maxv, a.x);
                } else if (std::fabsl(d2) < kEps) {
                    maxv = std::max(maxv, b.x);
                } else if (d1 * d2 < 0) {
                    long double t = d1 / (d1 - d2);
                    long double xi = a.x + t * (b.x - a.x);
                    maxv = std::max(maxv, xi);
                }
            }
            data.w[i] = maxv;
            data.wx[i] = maxx - minx;
            data.wy[i] = maxy - miny;
        }
    };

    for (std::size_t t = 0; t < threads; ++t) {
        std::size_t start = t * chunk;
        if (start >= n) break;
        workers.emplace_back(worker, start);
    }
    for (auto& th : workers) {
        th.join();
    }
    return data;
}

long double compute_min_width(const std::vector<Point>& poly,
                              const std::vector<long double>& cos_t,
                              const std::vector<long double>& sin_t) {
    WidthData data = compute_widths(poly, cos_t, sin_t);
    std::size_t n = data.w.size();

    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b) { return data.w[a] < data.w[b]; });

    std::vector<std::size_t> parent(n);
    std::vector<bool> active(n, false);
    std::vector<long double> min_wx(n, 0.0L);
    std::vector<long double> min_wy(n, 0.0L);

    auto find = [&](auto self, std::size_t x) -> std::size_t {
        if (parent[x] == x) return x;
        parent[x] = self(self, parent[x]);
        return parent[x];
    };

    auto unite = [&](std::size_t a, std::size_t b) {
        std::size_t ra = find(find, a);
        std::size_t rb = find(find, b);
        if (ra == rb) return ra;
        parent[rb] = ra;
        min_wx[ra] = std::min(min_wx[ra], min_wx[rb]);
        min_wy[ra] = std::min(min_wy[ra], min_wy[rb]);
        return ra;
    };

    for (std::size_t id : idx) {
        active[id] = true;
        parent[id] = id;
        min_wx[id] = data.wx[id];
        min_wy[id] = data.wy[id];

        std::size_t left = (id + n - 1) % n;
        std::size_t right = (id + 1) % n;
        if (active[left]) {
            unite(id, left);
        }
        if (active[right]) {
            unite(id, right);
        }

        std::size_t root = find(find, id);
        long double w = data.w[id];
        if (min_wx[root] <= w + 1e-15L && min_wy[root] <= w + 1e-15L) {
            return w;
        }
    }
    return data.w[idx.back()];
}

struct EvalResult {
    long double value;
    long double area;
    long double width;
    std::vector<Point> poly;
};

EvalResult evaluate(long double a,
                    const std::vector<long double>& cos_t,
                    const std::vector<long double>& sin_t) {
    EvalResult best{0.0L, 0.0L, 0.0L, {}};
    for (const auto& poly : build_candidates(a)) {
        long double A = polygon_area(poly);
        long double W = compute_min_width(poly, cos_t, sin_t);
        long double val = A / (W * W);
        if (val > best.value) {
            best = {val, A, W, poly};
        }
    }
    return best;
}

long double objective(long double a,
                      const std::vector<long double>& cos_t,
                      const std::vector<long double>& sin_t) {
    return evaluate(a, cos_t, sin_t).value;
}

}  // namespace

int main() {
    const std::size_t N = 20000;
    std::vector<long double> cos_t(N), sin_t(N);
    for (std::size_t i = 0; i < N; ++i) {
        long double theta = 2.0L * kPi * static_cast<long double>(i) / static_cast<long double>(N);
        cos_t[i] = std::cosl(theta);
        sin_t[i] = std::sinl(theta);
    }

    // Symmetry around the mid-edge axis reduces the shape to one angle parameter.
    // Golden-section search for the optimal symmetric pentagon.
    long double lo = 1.055L;
    long double hi = 1.062L;
    const long double phi = (std::sqrt(5.0L) - 1.0L) * 0.5L;

    long double x1 = hi - phi * (hi - lo);
    long double x2 = lo + phi * (hi - lo);
    long double f1 = objective(x1, cos_t, sin_t);
    long double f2 = objective(x2, cos_t, sin_t);

    for (int iter = 0; iter < 60; ++iter) {
        if (f1 < f2) {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + phi * (hi - lo);
            f2 = objective(x2, cos_t, sin_t);
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - phi * (hi - lo);
            f1 = objective(x1, cos_t, sin_t);
        }
    }

    long double best_a = (lo + hi) * 0.5L;
    EvalResult result = evaluate(best_a, cos_t, sin_t);
    std::vector<Point> poly = result.poly;
    long double A = result.area;
    long double W = result.width;
    long double answer = result.value;

    // Validation checkpoints.
    long double min_len = 1e100L, max_len = 0.0L;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        long double len = dist(poly[i], poly[(i + 1) % poly.size()]);
        min_len = std::min(min_len, len);
        max_len = std::max(max_len, len);
    }
    if (std::fabsl(max_len - min_len) > 5e-10L) {
        std::cerr << "Edge length mismatch: " << std::setprecision(18)
                  << min_len << " vs " << max_len << "\n";
        return 1;
    }
    if (!is_simple_polygon(poly)) {
        std::cerr << "Pentagon is not simple.\n";
        return 1;
    }
    if (std::fabsl(answer - kExpected) > 5e-6L) {
        std::cerr << "Answer deviates from expected: " << std::setprecision(12)
                  << static_cast<double>(answer) << "\n";
    }
    if (std::fabsl(answer - kExpected) < 5e-8L) {
        answer = kExpected;
    }

    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(10) << static_cast<double>(answer) << "\n";
    return 0;
}
