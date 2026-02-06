#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr int kTargetPointCount = 500;
constexpr i64 kMod = 50'515'093;

struct Point {
    int x;
    int y;

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};

inline i64 cross(const Point& a, const Point& b, const Point& c) {
    return static_cast<i64>(b.x - a.x) * (c.y - a.y) -
           static_cast<i64>(b.y - a.y) * (c.x - a.x);
}

inline i64 dist2(const Point& a, const Point& b) {
    const i64 dx = static_cast<i64>(a.x) - b.x;
    const i64 dy = static_cast<i64>(a.y) - b.y;
    return dx * dx + dy * dy;
}

std::vector<Point> generate_points(int count) {
    std::vector<int> t_values;
    t_values.reserve(static_cast<std::size_t>(count) * 2);

    i64 s = 290'797;
    for (int i = 0; i < count * 2; ++i) {
        s = (s * s) % kMod;
        t_values.push_back(static_cast<int>(s % 2000) - 1000);
    }

    std::vector<Point> points;
    points.reserve(count);
    for (int i = 0; i < count; ++i) {
        points.push_back(Point{t_values[2 * i], t_values[2 * i + 1]});
    }

    return points;
}

bool point_on_segment_inclusive(const Point& a, const Point& b, const Point& p) {
    if (cross(a, b, p) != 0) {
        return false;
    }
    const int min_x = std::min(a.x, b.x);
    const int max_x = std::max(a.x, b.x);
    const int min_y = std::min(a.y, b.y);
    const int max_y = std::max(a.y, b.y);
    return p.x >= min_x && p.x <= max_x && p.y >= min_y && p.y <= max_y;
}

struct GeometryCache {
    int n = 0;
    int word_count = 0;
    std::vector<u64> left_bits;   // left_bits[(a*n + b)*word_count + w]
    std::vector<std::uint8_t> has_inner_collinear; // has inner point on open segment a-b
};

GeometryCache precompute_geometry(const std::vector<Point>& points, unsigned threads) {
    GeometryCache cache;
    cache.n = static_cast<int>(points.size());
    cache.word_count = (cache.n + 63) / 64;
    cache.left_bits.assign(static_cast<std::size_t>(cache.n) * cache.n * cache.word_count, 0);
    cache.has_inner_collinear.assign(static_cast<std::size_t>(cache.n) * cache.n, 0);

    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }
    threads = std::max<unsigned>(1, std::min<unsigned>(threads, static_cast<unsigned>(cache.n)));

    auto worker = [&](int a_begin, int a_end) {
        for (int a = a_begin; a < a_end; ++a) {
            for (int b = 0; b < cache.n; ++b) {
                if (a == b) {
                    continue;
                }

                u64* bits = &cache.left_bits[(static_cast<std::size_t>(a) * cache.n + b) *
                                             cache.word_count];
                bool has_inner = false;

                for (int p = 0; p < cache.n; ++p) {
                    if (p == a || p == b) {
                        continue;
                    }

                    const i64 cr = cross(points[a], points[b], points[p]);
                    if (cr > 0) {
                        bits[p >> 6] |= (1ULL << (p & 63));
                    } else if (!has_inner && cr == 0 &&
                               point_on_segment_inclusive(points[a], points[b], points[p])) {
                        has_inner = true;
                    }
                }

                cache.has_inner_collinear[static_cast<std::size_t>(a) * cache.n + b] =
                    static_cast<std::uint8_t>(has_inner);
            }
        }
    };

    if (threads == 1 || cache.n < 64) {
        worker(0, cache.n);
        return cache;
    }

    std::vector<std::thread> pool;
    pool.reserve(threads);

    const int chunk = (cache.n + static_cast<int>(threads) - 1) / static_cast<int>(threads);
    for (unsigned tid = 0; tid < threads; ++tid) {
        const int begin = static_cast<int>(tid) * chunk;
        const int end = std::min(cache.n, begin + chunk);
        if (begin >= end) {
            continue;
        }
        pool.emplace_back(worker, begin, end);
    }

    for (std::thread& th : pool) {
        th.join();
    }

    return cache;
}

inline bool triangle_is_empty_strict(int a,
                                     int b,
                                     int c,
                                     const GeometryCache& cache) {
    const u64* ab = &cache.left_bits[(static_cast<std::size_t>(a) * cache.n + b) *
                                     cache.word_count];
    const u64* bc = &cache.left_bits[(static_cast<std::size_t>(b) * cache.n + c) *
                                     cache.word_count];
    const u64* ca = &cache.left_bits[(static_cast<std::size_t>(c) * cache.n + a) *
                                     cache.word_count];

    for (int w = 0; w < cache.word_count; ++w) {
        if ((ab[w] & bc[w] & ca[w]) != 0) {
            return false;
        }
    }
    return true;
}

i64 solve_for_anchor(int anchor,
                     const std::vector<Point>& points,
                     const GeometryCache& cache) {
    const int n = static_cast<int>(points.size());

    std::vector<int> verts;
    verts.reserve(n - 1);
    for (int idx = 0; idx < n; ++idx) {
        if (idx == anchor) {
            continue;
        }
        if (points[idx].y > points[anchor].y ||
            (points[idx].y == points[anchor].y && points[idx].x > points[anchor].x)) {
            verts.push_back(idx);
        }
    }

    const int m = static_cast<int>(verts.size());
    if (m < 2) {
        return 0;
    }

    auto angle_cmp = [&](int lhs, int rhs) {
        const i64 cr = cross(points[anchor], points[lhs], points[rhs]);
        if (cr != 0) {
            return cr > 0;
        }
        const i64 d_lhs = dist2(points[anchor], points[lhs]);
        const i64 d_rhs = dist2(points[anchor], points[rhs]);
        if (d_lhs != d_rhs) {
            return d_lhs < d_rhs;
        }
        return lhs < rhs;
    };
    std::sort(verts.begin(), verts.end(), angle_cmp);

    std::vector<i64> dp(static_cast<std::size_t>(m) * m, -1);
    std::vector<std::vector<int>> incoming(m);

    i64 best = 0;

    for (int j = 0; j < m - 1; ++j) {
        const int vj = verts[j];
        const bool can_extend_over_j =
            !cache.has_inner_collinear[static_cast<std::size_t>(anchor) * n + vj];

        for (int k = j + 1; k < m; ++k) {
            const int vk = verts[k];
            const i64 tri_area2 = cross(points[anchor], points[vj], points[vk]);
            if (tri_area2 <= 0) {
                continue;
            }

            if (!triangle_is_empty_strict(anchor, vj, vk, cache)) {
                continue;
            }

            i64 current = tri_area2; // base triangle (anchor, vj, vk)

            if (can_extend_over_j) {
                for (const int h : incoming[j]) {
                    const int vh = verts[h];
                    if (cross(points[vh], points[vj], points[vk]) <= 0) {
                        continue;
                    }
                    const i64 candidate = dp[static_cast<std::size_t>(h) * m + j] + tri_area2;
                    if (candidate > current) {
                        current = candidate;
                    }
                }
            }

            dp[static_cast<std::size_t>(j) * m + k] = current;
            incoming[k].push_back(j);
            if (current > best) {
                best = current;
            }
        }
    }

    return best;
}

i64 solve_max_area2(const std::vector<Point>& points,
                    bool allow_multithreading,
                    unsigned requested_threads = 0) {
    if (points.size() < 3) {
        return 0;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    if (!allow_multithreading || points.size() < 120) {
        threads = 1;
    }

    const GeometryCache cache = precompute_geometry(points, threads);

    if (threads == 1) {
        i64 best = 0;
        for (int anchor = 0; anchor < static_cast<int>(points.size()); ++anchor) {
            const i64 local = solve_for_anchor(anchor, points, cache);
            if (local > best) {
                best = local;
            }
        }
        return best;
    }

    std::atomic<int> next_anchor{0};
    std::vector<i64> partial_best(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            i64 local_best = 0;
            while (true) {
                const int anchor = next_anchor.fetch_add(1, std::memory_order_relaxed);
                if (anchor >= static_cast<int>(points.size())) {
                    break;
                }
                const i64 value = solve_for_anchor(anchor, points, cache);
                if (value > local_best) {
                    local_best = value;
                }
            }
            partial_best[tid] = local_best;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    i64 best = 0;
    for (const i64 value : partial_best) {
        if (value > best) {
            best = value;
        }
    }
    return best;
}

std::vector<int> convex_hull_indices(std::vector<int> ids, const std::vector<Point>& points) {
    std::sort(ids.begin(), ids.end(), [&](int lhs, int rhs) {
        if (points[lhs].x != points[rhs].x) {
            return points[lhs].x < points[rhs].x;
        }
        if (points[lhs].y != points[rhs].y) {
            return points[lhs].y < points[rhs].y;
        }
        return lhs < rhs;
    });

    std::vector<int> lower;
    for (const int idx : ids) {
        while (lower.size() >= 2 &&
               cross(points[lower[lower.size() - 2]], points[lower.back()], points[idx]) <= 0) {
            lower.pop_back();
        }
        lower.push_back(idx);
    }

    std::vector<int> upper;
    for (int i = static_cast<int>(ids.size()) - 1; i >= 0; --i) {
        const int idx = ids[i];
        while (upper.size() >= 2 &&
               cross(points[upper[upper.size() - 2]], points[upper.back()], points[idx]) <= 0) {
            upper.pop_back();
        }
        upper.push_back(idx);
    }

    if (lower.size() <= 1 || upper.size() <= 1) {
        return {};
    }

    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());

    if (lower.size() < 3) {
        return {};
    }
    return lower;
}

i64 polygon_area2(const std::vector<int>& hull, const std::vector<Point>& points) {
    i64 twice_area = 0;
    for (std::size_t i = 0; i < hull.size(); ++i) {
        const Point& a = points[hull[i]];
        const Point& b = points[hull[(i + 1) % hull.size()]];
        twice_area += static_cast<i64>(a.x) * b.y - static_cast<i64>(a.y) * b.x;
    }
    if (twice_area < 0) {
        twice_area = -twice_area;
    }
    return twice_area;
}

bool point_strictly_inside_convex_polygon(const Point& p,
                                          const std::vector<int>& hull,
                                          const std::vector<Point>& points) {
    bool on_boundary = false;

    for (std::size_t i = 0; i < hull.size(); ++i) {
        const Point& a = points[hull[i]];
        const Point& b = points[hull[(i + 1) % hull.size()]];
        const i64 cr = cross(a, b, p);
        if (cr < 0) {
            return false;
        }
        if (cr == 0 && point_on_segment_inclusive(a, b, p)) {
            on_boundary = true;
        }
    }

    return !on_boundary;
}

i64 brute_force_max_area2(const std::vector<Point>& points) {
    const int n = static_cast<int>(points.size());
    if (n >= 20) {
        return -1; // safety guard: not intended for larger N.
    }

    i64 best = 0;
    const int total_masks = 1 << n;
    for (int mask = 0; mask < total_masks; ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) < 3) {
            continue;
        }

        std::vector<int> subset;
        subset.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {
                subset.push_back(i);
            }
        }

        const std::vector<int> hull = convex_hull_indices(subset, points);
        if (hull.size() != subset.size()) {
            continue;
        }

        const i64 area2 = polygon_area2(hull, points);
        if (area2 <= best) {
            continue;
        }

        bool ok = true;
        for (int p = 0; p < n; ++p) {
            if (mask & (1 << p)) {
                continue;
            }
            if (point_strictly_inside_convex_polygon(points[p], hull, points)) {
                ok = false;
                break;
            }
        }

        if (ok) {
            best = area2;
        }
    }

    return best;
}

bool run_validation_checkpoints() {
    {
        const std::vector<Point> p = generate_points(3);
        const std::vector<Point> expected = {
            {527, 144},
            {-488, 732},
            {-454, -947},
        };
        if (p != expected) {
            std::cerr << "Generator checkpoint failed.\n";
            return false;
        }
    }

    {
        const std::vector<Point> p = generate_points(10);
        const i64 fast = solve_max_area2(p, false, 1);
        const i64 brute = brute_force_max_area2(p);
        if (fast != brute) {
            std::cerr << "Brute-force checkpoint failed for N=10: got " << fast
                      << ", expected " << brute << "\n";
            return false;
        }
    }

    {
        const std::vector<Point> p = generate_points(20);
        constexpr i64 kExpectedArea2 = 2'099'389; // 1049694.5
        const i64 got = solve_max_area2(p, false, 1);
        if (got != kExpectedArea2) {
            std::cerr << "Sample checkpoint failed for N=20: got " << got
                      << ", expected " << kExpectedArea2 << "\n";
            return false;
        }
    }

    return true;
}

void print_area_with_single_decimal(i64 area2) {
    std::cout << (area2 / 2) << ((area2 & 1) ? ".5" : ".0") << '\n';
}

} // namespace

int main() {
    if (!run_validation_checkpoints()) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    const std::vector<Point> points = generate_points(kTargetPointCount);
    const i64 answer_area2 = solve_max_area2(points, true, threads);
    print_area_with_single_decimal(answer_area2);
    return 0;
}
