#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;

namespace {

constexpr long double kSqrt2 = 1.4142135623730950488016887242097L;

struct Num {
    long long a = 0;
    long long b = 0;  // value = a + b*sqrt(2)

    Num() = default;
    Num(long long a_, long long b_) : a(a_), b(b_) {}

    Num operator+(const Num& other) const { return {a + other.a, b + other.b}; }
    Num operator-(const Num& other) const { return {a - other.a, b - other.b}; }
    Num operator*(const Num& other) const {
        // (a + b*sqrt(2))(c + d*sqrt(2)) = (ac + 2bd) + (ad + bc)*sqrt(2)
        return {a * other.a + 2 * b * other.b, a * other.b + b * other.a};
    }

    bool operator==(const Num& other) const { return a == other.a && b == other.b; }
    bool is_zero() const { return a == 0 && b == 0; }

    long double value() const { return static_cast<long double>(a) + static_cast<long double>(b) * kSqrt2; }
};

struct Point {
    Num x;
    Num y;

    bool operator==(const Point& other) const { return x == other.x && y == other.y; }
};

struct PointKey {
    long long xa = 0;
    long long xb = 0;
    long long ya = 0;
    long long yb = 0;

    bool operator==(const PointKey& other) const {
        return xa == other.xa && xb == other.xb && ya == other.ya && yb == other.yb;
    }
};

PointKey key_of(const Point& p) { return {p.x.a, p.x.b, p.y.a, p.y.b}; }

Point point_from_key(const PointKey& k) { return {{k.xa, k.xb}, {k.ya, k.yb}}; }

struct PointKeyHash {
    size_t operator()(const PointKey& k) const {
        size_t h = std::hash<long long>{}(k.xa);
        h ^= std::hash<long long>{}(k.xb) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<long long>{}(k.ya) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<long long>{}(k.yb) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct EdgeKey {
    PointKey a;
    PointKey b;

    bool operator==(const EdgeKey& other) const { return a == other.a && b == other.b; }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const {
        size_t h = PointKeyHash{}(k.a);
        h ^= PointKeyHash{}(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct EdgeVal {
    Point from;
    Point to;
};

int sign_num(const Num& v) {
    if (v.is_zero()) return 0;
    return v.value() < 0 ? -1 : 1;
}

int cmp_num(const Num& a, const Num& b) { return sign_num(a - b); }

bool point_less_geom(const Point& a, const Point& b) {
    int cy = cmp_num(a.y, b.y);
    if (cy != 0) return cy < 0;
    return cmp_num(a.x, b.x) < 0;
}

bool pointkey_less(const PointKey& a, const PointKey& b) {
    if (a.xa != b.xa) return a.xa < b.xa;
    if (a.xb != b.xb) return a.xb < b.xb;
    if (a.ya != b.ya) return a.ya < b.ya;
    return a.yb < b.yb;
}

Point operator+(const Point& a, const Point& b) { return {a.x + b.x, a.y + b.y}; }
Point operator-(const Point& a, const Point& b) { return {a.x - b.x, a.y - b.y}; }

Num cross(const Point& a, const Point& b, const Point& c) {
    Point ab = b - a;
    Point ac = c - a;
    Num term1 = ab.x * ac.y;
    Num term2 = ab.y * ac.x;
    return term1 - term2;
}

Num min_num(const Num& a, const Num& b) { return cmp_num(a, b) <= 0 ? a : b; }
Num max_num(const Num& a, const Num& b) { return cmp_num(a, b) >= 0 ? a : b; }

bool on_segment(const Point& a, const Point& b, const Point& p) {
    if (!cross(a, b, p).is_zero()) return false;
    if (cmp_num(p.x, min_num(a.x, b.x)) < 0 || cmp_num(p.x, max_num(a.x, b.x)) > 0) return false;
    if (cmp_num(p.y, min_num(a.y, b.y)) < 0 || cmp_num(p.y, max_num(a.y, b.y)) > 0) return false;
    return true;
}

int orient(const Point& a, const Point& b, const Point& c) { return sign_num(cross(a, b, c)); }

bool seg_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    int o1 = orient(a, b, c);
    int o2 = orient(a, b, d);
    int o3 = orient(c, d, a);
    int o4 = orient(c, d, b);
    if (o1 == 0 && on_segment(a, b, c)) return true;
    if (o2 == 0 && on_segment(a, b, d)) return true;
    if (o3 == 0 && on_segment(c, d, a)) return true;
    if (o4 == 0 && on_segment(c, d, b)) return true;
    return (o1 != o2 && o3 != o4);
}

bool point_in_polygon(const Point& p, const vector<Point>& poly) {
    long double px = p.x.value();
    long double py = p.y.value();
    bool inside = false;
    size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        const Point& a = poly[i];
        const Point& b = poly[(i + 1) % n];
        if (on_segment(a, b, p)) return true;
        long double ay = a.y.value();
        long double by = b.y.value();
        if ((ay > py) != (by > py)) {
            long double ax = a.x.value();
            long double bx = b.x.value();
            long double xin = ax + (bx - ax) * (py - ay) / (by - ay);
            if (xin > px) inside = !inside;
        }
    }
    return inside;
}

EdgeKey make_edge_key(const Point& a, const Point& b) {
    PointKey ka = key_of(a);
    PointKey kb = key_of(b);
    if (pointkey_less(ka, kb)) return {ka, kb};
    return {kb, ka};
}

vector<Point> build_vertices(const vector<uint8_t>& steps, const array<Point, 8>& dirs) {
    vector<Point> pts;
    pts.reserve(steps.size() + 1);
    Point cur{{0, 0}, {0, 0}};
    pts.push_back(cur);
    for (uint8_t d : steps) {
        cur = cur + dirs[d];
        pts.push_back(cur);
    }
    return pts;
}

long double signed_area(const vector<uint8_t>& steps, const array<Point, 8>& dirs) {
    if (steps.empty()) return 0.0L;
    vector<Point> pts = build_vertices(steps, dirs);
    // drop repeated last point
    pts.pop_back();
    long double area = 0.0L;
    for (size_t i = 0; i < pts.size(); ++i) {
        const Point& a = pts[i];
        const Point& b = pts[(i + 1) % pts.size()];
        area += a.x.value() * b.y.value() - b.x.value() * a.y.value();
    }
    return 0.5L * area;
}

vector<uint8_t> to_ccw(const vector<uint8_t>& steps) {
    vector<uint8_t> out;
    out.reserve(steps.size());
    for (size_t i = steps.size(); i > 0; --i) {
        out.push_back(static_cast<uint8_t>((steps[i - 1] + 4) & 7));
    }
    return out;
}

unordered_map<EdgeKey, EdgeVal, EdgeKeyHash> edges_from_steps(
    const vector<uint8_t>& steps, const array<Point, 8>& dirs) {
    unordered_map<EdgeKey, EdgeVal, EdgeKeyHash> edges;
    edges.reserve(steps.size() * 2);
    Point cur{{0, 0}, {0, 0}};
    for (uint8_t d : steps) {
        Point nxt = cur + dirs[d];
        edges.emplace(make_edge_key(cur, nxt), EdgeVal{cur, nxt});
        cur = nxt;
    }
    return edges;
}

bool extract_cycles_from_edges(const unordered_map<EdgeKey, EdgeVal, EdgeKeyHash>& edge_map,
                               const array<Point, 8>& dirs,
                               vector<vector<uint8_t>>& cycles) {
    // Walk each directed boundary cycle by keeping the interior on the left
    // (smallest CCW turn from the incoming direction), even when components touch.
    cycles.clear();
    if (edge_map.empty()) return true;

    struct EdgeInfo {
        PointKey from;
        PointKey to;
        int dir = -1;
        bool used = false;
    };

    vector<EdgeInfo> edges;
    edges.reserve(edge_map.size());
    unordered_map<PointKey, vector<int>, PointKeyHash> outgoing;
    outgoing.reserve(edge_map.size() * 2);

    for (const auto& kv : edge_map) {
        const EdgeVal& e = kv.second;
        Point diff = e.to - e.from;
        int dir_idx = -1;
        for (int i = 0; i < 8; ++i) {
            if (diff.x == dirs[i].x && diff.y == dirs[i].y) {
                dir_idx = i;
                break;
            }
        }
        if (dir_idx < 0) return false;
        int id = static_cast<int>(edges.size());
        edges.push_back({key_of(e.from), key_of(e.to), dir_idx, false});
        outgoing[edges.back().from].push_back(id);
    }

    auto pick_start = [&]() -> int {
        bool init = false;
        Point min_point;
        int best_id = -1;
        for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
            if (edges[i].used) continue;
            Point p = point_from_key(edges[i].from);
            if (!init || point_less_geom(p, min_point)) {
                min_point = p;
                best_id = i;
                init = true;
            }
        }
        return best_id;
    };

    size_t edges_used = 0;
    while (edges_used < edges.size()) {
        int start_id = pick_start();
        if (start_id < 0) return false;
        int cur_id = start_id;
        vector<uint8_t> steps;
        while (true) {
            EdgeInfo& cur = edges[cur_id];
            if (cur.used) return false;
            cur.used = true;
            ++edges_used;
            steps.push_back(static_cast<uint8_t>(cur.dir));

            PointKey head = cur.to;
            auto it_out = outgoing.find(head);
            if (it_out == outgoing.end()) return false;

            int best_next = -1;
            int best_delta = 9;
            for (int cand_id : it_out->second) {
                if (edges[cand_id].used && cand_id != start_id) continue;
                int delta = (edges[cand_id].dir - cur.dir + 8) & 7;
                if (delta < best_delta) {
                    best_delta = delta;
                    best_next = cand_id;
                }
            }
            if (best_next < 0) return false;
            if (best_next == start_id) {
                break;
            }
            cur_id = best_next;
        }
        cycles.push_back(std::move(steps));
    }

    return edges_used == edges.size();
}

string steps_key(const vector<uint8_t>& steps) {
    string key;
    key.resize(steps.size());
    for (size_t i = 0; i < steps.size(); ++i) key[i] = static_cast<char>(steps[i]);
    return key;
}

struct Solver {
    array<Point, 8> dirs;
    unordered_map<string, cpp_int> memo;

    Solver() {
        // Directions are scaled by 2 to avoid halves: (x, y) = (a + b*sqrt(2), c + d*sqrt(2))
        dirs[0] = {{2, 0}, {0, 0}};   // E
        dirs[1] = {{0, 1}, {0, 1}};   // NE
        dirs[2] = {{0, 0}, {2, 0}};   // N
        dirs[3] = {{0, -1}, {0, 1}};  // NW
        dirs[4] = {{-2, 0}, {0, 0}};  // W
        dirs[5] = {{0, -1}, {0, -1}}; // SW
        dirs[6] = {{0, 0}, {-2, 0}};  // S
        dirs[7] = {{0, 1}, {0, -1}};  // SE
    }

    cpp_int dfs(const string& key) {
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;
        if (key.empty()) return cpp_int(1);

        vector<uint8_t> steps;
        steps.reserve(key.size());
        for (unsigned char c : key) steps.push_back(static_cast<uint8_t>(c));

        vector<Point> pts = build_vertices(steps, dirs);
        size_t n = steps.size();

        // Find the lowest (y, then x) vertex.
        size_t min_idx = 0;
        Point min_point = pts[0];
        for (size_t i = 0; i < n; ++i) {
            if (point_less_geom(pts[i], min_point)) {
                min_point = pts[i];
                min_idx = i;
            }
        }

        const Point p = min_point;
        const uint8_t d = steps[min_idx];
        vector<pair<Point, Point>> boundary_edges;
        boundary_edges.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            boundary_edges.emplace_back(pts[i], pts[i + 1]);
        }

        cpp_int total = 0;
        auto base_edges = edges_from_steps(steps, dirs);

        for (int delta = 1; delta <= 3; ++delta) {
            uint8_t vd = static_cast<uint8_t>((d + delta) & 7);
            Point u = dirs[d];
            Point v = dirs[vd];
            Point p1 = p;
            Point p2 = p + u;
            Point p3 = p2 + v;
            Point p4 = p + v;

            if (!point_in_polygon(p4, pts) || !point_in_polygon(p3, pts)) continue;

            bool ok = true;
            array<pair<Point, Point>, 4> tile_edges = {{
                {p1, p2}, {p2, p3}, {p3, p4}, {p4, p1},
            }};
            for (const auto& te : tile_edges) {
                for (const auto& be : boundary_edges) {
                    const Point& a = te.first;
                    const Point& b = te.second;
                    const Point& c = be.first;
                    const Point& d2 = be.second;
                    if ((a == c && b == d2) || (a == d2 && b == c)) continue;
                    if (a == c || a == d2 || b == c || b == d2) continue;
                    if (seg_intersect(a, b, c, d2)) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) break;
            }
            if (!ok) continue;

            auto edges = base_edges;
            array<pair<Point, Point>, 4> cw_edges = {{
                {p1, p4}, {p4, p3}, {p3, p2}, {p2, p1},
            }};
            for (const auto& e : cw_edges) {
                EdgeKey ek = make_edge_key(e.first, e.second);
                auto it_edge = edges.find(ek);
                if (it_edge != edges.end()) {
                    edges.erase(it_edge);
                } else {
                    edges.emplace(ek, EdgeVal{e.first, e.second});
                }
            }

            vector<vector<uint8_t>> cycles;
            if (!extract_cycles_from_edges(edges, dirs, cycles)) continue;
            if (cycles.empty()) {
                total += 1;
                continue;
            }
            cpp_int ways = 1;
            for (auto& cycle : cycles) {
                if (signed_area(cycle, dirs) < 0) cycle = to_ccw(cycle);
                ways *= dfs(steps_key(cycle));
            }
            total += ways;
        }

        memo.emplace(key, total);
        return total;
    }

    cpp_int solve(int a, int b) {
        vector<uint8_t> steps;
        steps.reserve(8 * (a + b));
        array<int, 8> lens = {a, b, a, b, a, b, a, b};
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < lens[i]; ++j) steps.push_back(static_cast<uint8_t>(i));
        }
        return dfs(steps_key(steps));
    }
};

bool run_validation(unsigned threads) {
    struct Case {
        int a;
        int b;
        const char* expected;
    };
    const Case cases[] = {
        {1, 1, "8"},
        {2, 1, "76"},
        {3, 2, "456572"},
    };

    bool ok = true;
    for (const auto& c : cases) {
        Solver solver;
        cpp_int got = solver.solve(c.a, c.b);
        std::string got_str = got.convert_to<std::string>();
        if (got_str != c.expected) {
            cerr << "Validation failed for O_" << c.a << "," << c.b
                 << ": got " << got_str << ", expected " << c.expected << "\n";
            ok = false;
        }
    }
    if (ok) cerr << "Validation checkpoints passed.\n";
    (void)threads;
    return ok;
}

cpp_int solve_parallel(int a, int b, unsigned threads) {
    Solver base;
    vector<uint8_t> steps;
    steps.reserve(8 * (a + b));
    array<int, 8> lens = {a, b, a, b, a, b, a, b};
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < lens[i]; ++j) steps.push_back(static_cast<uint8_t>(i));
    }

    // Build initial state and split on the first edge placements.
    vector<Point> pts = build_vertices(steps, base.dirs);
    size_t n = steps.size();
    size_t min_idx = 0;
    Point min_point = pts[0];
    for (size_t i = 0; i < n; ++i) {
        if (point_less_geom(pts[i], min_point)) {
            min_point = pts[i];
            min_idx = i;
        }
    }
    Point p = min_point;
    uint8_t d = steps[min_idx];
    vector<pair<Point, Point>> boundary_edges;
    boundary_edges.reserve(n);
    for (size_t i = 0; i < n; ++i) boundary_edges.emplace_back(pts[i], pts[i + 1]);

    auto base_edges = edges_from_steps(steps, base.dirs);
    vector<string> sub_tasks;
    vector<cpp_int> fixed_tasks;
    for (int delta = 1; delta <= 3; ++delta) {
        uint8_t vd = static_cast<uint8_t>((d + delta) & 7);
        Point u = base.dirs[d];
        Point v = base.dirs[vd];
        Point p1 = p;
        Point p2 = p + u;
        Point p3 = p2 + v;
        Point p4 = p + v;
        if (!point_in_polygon(p4, pts) || !point_in_polygon(p3, pts)) continue;

        bool ok = true;
        array<pair<Point, Point>, 4> tile_edges = {{
            {p1, p2}, {p2, p3}, {p3, p4}, {p4, p1},
        }};
        for (const auto& te : tile_edges) {
            for (const auto& be : boundary_edges) {
                const Point& a0 = te.first;
                const Point& b0 = te.second;
                const Point& c0 = be.first;
                const Point& d0 = be.second;
                if ((a0 == c0 && b0 == d0) || (a0 == d0 && b0 == c0)) continue;
                if (a0 == c0 || a0 == d0 || b0 == c0 || b0 == d0) continue;
                if (seg_intersect(a0, b0, c0, d0)) {
                    ok = false;
                    break;
                }
            }
            if (!ok) break;
        }
        if (!ok) continue;

        auto edges = base_edges;
        array<pair<Point, Point>, 4> cw_edges = {{
            {p1, p4}, {p4, p3}, {p3, p2}, {p2, p1},
        }};
        for (const auto& e : cw_edges) {
            EdgeKey ek = make_edge_key(e.first, e.second);
            auto it_edge = edges.find(ek);
            if (it_edge != edges.end()) {
                edges.erase(it_edge);
            } else {
                edges.emplace(ek, EdgeVal{e.first, e.second});
            }
        }
        vector<vector<uint8_t>> cycles;
        if (!extract_cycles_from_edges(edges, base.dirs, cycles)) continue;
        if (cycles.empty()) {
            fixed_tasks.push_back(1);
            continue;
        }
        if (cycles.size() == 1) {
            if (signed_area(cycles[0], base.dirs) < 0) cycles[0] = to_ccw(cycles[0]);
            sub_tasks.push_back(steps_key(cycles[0]));
            continue;
        }
        // Multiple components: solve independently and multiply here.
        cpp_int ways = 1;
        for (auto& cycle : cycles) {
            if (signed_area(cycle, base.dirs) < 0) cycle = to_ccw(cycle);
            Solver solver;
            ways *= solver.dfs(steps_key(cycle));
        }
        fixed_tasks.push_back(ways);
    }

    if (threads <= 1 || sub_tasks.size() <= 1) {
        Solver solver;
        cpp_int total = 0;
        for (const auto& val : fixed_tasks) total += val;
        for (const auto& key : sub_tasks) total += solver.dfs(key);
        return total;
    }

    vector<future<cpp_int>> futures;
    futures.reserve(sub_tasks.size());
    for (const auto& key : sub_tasks) {
        futures.push_back(std::async(std::launch::async, [key]() {
            Solver solver;
            return solver.dfs(key);
        }));
    }

    cpp_int total = 0;
    for (const auto& val : fixed_tasks) total += val;
    for (auto& f : futures) total += f.get();
    return total;
}

} // namespace

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int a = 4;
    int b = 2;
    unsigned threads = thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = min(threads, 4u);
    bool validate = true;

    // Optional CLI: ./a.out [a] [b] [threads] [validate(0/1)]
    if (argc >= 2) a = stoi(argv[1]);
    if (argc >= 3) b = stoi(argv[2]);
    if (argc >= 4) threads = static_cast<unsigned>(stoul(argv[3]));
    if (argc >= 5) validate = (stoi(argv[4]) != 0);

    if (validate && !run_validation(min(threads, 2u))) return 1;

    cpp_int answer = solve_parallel(a, b, threads);
    cout << answer << "\n";
    return 0;
}
