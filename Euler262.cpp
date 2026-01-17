#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct Point {
    double x = 0.0;
    double y = 0.0;
};

double height(double x, double y) {
    const double xx = x * x;
    const double yy = y * y;
    const double xy = x * y;
    const double base = 5000.0 - 0.005 * (xx + yy + xy) + 12.5 * (x + y);
    const double expo = -std::abs(0.000001 * (xx + yy) - 0.0015 * (x + y) + 0.7);
    return base * std::exp(expo);
}

double dist(const Point &a, const Point &b) {
    return std::hypot(a.x - b.x, a.y - b.y);
}

struct Key {
    int64_t x = 0;
    int64_t y = 0;

    bool operator==(const Key &other) const { return x == other.x && y == other.y; }
};

struct KeyHash {
    size_t operator()(const Key &k) const {
        uint64_t x = static_cast<uint64_t>(k.x);
        uint64_t y = static_cast<uint64_t>(k.y);
        return static_cast<size_t>((x * 0x9e3779b97f4a7c15ULL) ^ (y + (x << 6) + (x >> 2)));
    }
};

Key quantize(const Point &p, double scale) {
    return {static_cast<int64_t>(std::llround(p.x * scale)),
            static_cast<int64_t>(std::llround(p.y * scale))};
}

Point from_key(const Key &k, double scale) {
    return {static_cast<double>(k.x) / scale, static_cast<double>(k.y) / scale};
}

Point interp(const Point &p1, const Point &p2, double v1, double v2, double f) {
    double t;
    const double denom = v2 - v1;
    if (std::abs(denom) < 1e-14) {
        t = 0.5;
    } else {
        t = (f - v1) / denom;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
    }
    return {p1.x + t * (p2.x - p1.x), p1.y + t * (p2.y - p1.y)};
}

void compute_heights(std::vector<double> &heights, int n, double step, int threads) {
    auto worker = [&](int row_start, int row_end) {
        for (int i = row_start; i < row_end; ++i) {
            const double x = i * step;
            const int base = i * n;
            for (int j = 0; j < n; ++j) {
                const double y = j * step;
                heights[base + j] = height(x, y);
            }
        }
    };

    if (threads <= 1) {
        worker(0, n);
        return;
    }

    int use_threads = std::min(threads, n);
    int chunk = (n + use_threads - 1) / use_threads;
    std::vector<std::thread> pool;
    pool.reserve(use_threads);
    for (int t = 0; t < use_threads; ++t) {
        int start = t * chunk;
        int end = std::min(n, start + chunk);
        if (start >= end) break;
        pool.emplace_back(worker, start, end);
    }
    for (auto &th : pool) th.join();
}

double minimax_height(const std::vector<double> &heights, int n, int start_idx, int goal_idx) {
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> best(heights.size(), inf);

    using Node = std::pair<double, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    best[start_idx] = heights[start_idx];
    pq.emplace(best[start_idx], start_idx);

    const int dirs[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    while (!pq.empty()) {
        auto [cost, idx] = pq.top();
        pq.pop();
        if (cost != best[idx]) continue;
        if (idx == goal_idx) return cost;

        int i = idx / n;
        int j = idx - i * n;
        for (auto &d : dirs) {
            int ni = i + d[0];
            int nj = j + d[1];
            if (ni < 0 || ni >= n || nj < 0 || nj >= n) continue;
            int nidx = ni * n + nj;
            double next = std::max(cost, heights[nidx]);
            if (next < best[nidx]) {
                best[nidx] = next;
                pq.emplace(next, nidx);
            }
        }
    }
    return inf;
}

bool visible(const Point &a, const Point &b, double f, double step) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double dist_ab = std::hypot(dx, dy);
    if (dist_ab == 0.0) return true;
    const int samples = std::max(1, static_cast<int>(dist_ab / step));
    for (int s = 1; s < samples; ++s) {
        const double t = static_cast<double>(s) / samples;
        const double x = a.x + dx * t;
        const double y = a.y + dy * t;
        if (height(x, y) > f + 1e-10) return false;
    }
    return true;
}

struct EdgePair {
    int a = -1;
    int b = -1;
};

}  // namespace

int main(int argc, char **argv) {
    constexpr int kMaxCoord = 1600;
    double grid_step = 1.0;
    double vis_step = 0.5;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;

    bool verbose = false;
    bool validate = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            threads = std::max(1, std::atoi(argv[++i]));
        } else if (arg == "--grid-step" && i + 1 < argc) {
            grid_step = std::max(0.25, std::atof(argv[++i]));
        } else if (arg == "--vis-step" && i + 1 < argc) {
            vis_step = std::max(0.05, std::atof(argv[++i]));
        } else if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--validate") {
            validate = true;
        }
    }

    const int n = static_cast<int>(kMaxCoord / grid_step) + 1;
    std::vector<double> heights(static_cast<size_t>(n) * n, 0.0);
    compute_heights(heights, n, grid_step, threads);

    const Point A{200.0, 200.0};
    const Point B{1400.0, 1400.0};
    const int start_idx = static_cast<int>(A.x / grid_step) * n + static_cast<int>(A.y / grid_step);
    const int goal_idx = static_cast<int>(B.x / grid_step) * n + static_cast<int>(B.y / grid_step);

    const double f_min = minimax_height(heights, n, start_idx, goal_idx);
    if (!std::isfinite(f_min)) {
        std::cerr << "Failed to find a connecting elevation.\n";
        return 1;
    }

    if (f_min + 1e-9 < height(A.x, A.y) || f_min + 1e-9 < height(B.x, B.y)) {
        std::cerr << "Validation failed: f_min below endpoint height.\n";
        return 1;
    }

    if (verbose) {
        std::cerr << std::fixed << std::setprecision(6) << "f_min=" << f_min << "\n";
    }

    std::vector<uint8_t> inside(static_cast<size_t>(n) * n, 0);
    for (size_t idx = 0; idx < heights.size(); ++idx) {
        inside[idx] = heights[idx] > f_min ? 1 : 0;
    }

    static const int kCaseCount[16] = {
        0, 1, 1, 1, 1, 2, 1, 1,
        1, 1, 2, 1, 1, 1, 1, 0
    };
    static const EdgePair kCaseEdges[16][2] = {
        {{-1, -1}, {-1, -1}},  // 0
        {{3, 0}, {-1, -1}},    // 1
        {{0, 1}, {-1, -1}},    // 2
        {{3, 1}, {-1, -1}},    // 3
        {{1, 2}, {-1, -1}},    // 4
        {{3, 2}, {0, 1}},      // 5
        {{0, 2}, {-1, -1}},    // 6
        {{3, 2}, {-1, -1}},    // 7
        {{2, 3}, {-1, -1}},    // 8
        {{0, 2}, {-1, -1}},    // 9
        {{0, 1}, {2, 3}},      // 10
        {{1, 2}, {-1, -1}},    // 11
        {{3, 1}, {-1, -1}},    // 12
        {{0, 1}, {-1, -1}},    // 13
        {{3, 0}, {-1, -1}},    // 14
        {{-1, -1}, {-1, -1}}   // 15
    };

    std::vector<std::pair<Point, Point>> segments;
    segments.reserve(9000);

    for (int i = 0; i < n - 1; ++i) {
        const double x = i * grid_step;
        for (int j = 0; j < n - 1; ++j) {
            const double y = j * grid_step;

            const int idx_ll = i * n + j;
            const int idx_lr = (i + 1) * n + j;
            const int idx_ur = (i + 1) * n + (j + 1);
            const int idx_ul = i * n + (j + 1);

            const bool c0 = inside[idx_ll];
            const bool c1 = inside[idx_lr];
            const bool c2 = inside[idx_ur];
            const bool c3 = inside[idx_ul];
            const int mask = (c0 ? 1 : 0) | (c1 ? 2 : 0) | (c2 ? 4 : 0) | (c3 ? 8 : 0);
            if (mask == 0 || mask == 15) continue;

            const Point ll{x, y};
            const Point lr{x + grid_step, y};
            const Point ur{x + grid_step, y + grid_step};
            const Point ul{x, y + grid_step};

            const double v_ll = heights[idx_ll];
            const double v_lr = heights[idx_lr];
            const double v_ur = heights[idx_ur];
            const double v_ul = heights[idx_ul];

            Point edges[4];
            edges[0] = interp(ll, lr, v_ll, v_lr, f_min);
            edges[1] = interp(lr, ur, v_lr, v_ur, f_min);
            edges[2] = interp(ul, ur, v_ul, v_ur, f_min);
            edges[3] = interp(ll, ul, v_ll, v_ul, f_min);

            for (int s = 0; s < kCaseCount[mask]; ++s) {
                int e0 = kCaseEdges[mask][s].a;
                int e1 = kCaseEdges[mask][s].b;
                if (e0 < 0 || e1 < 0) continue;
                segments.emplace_back(edges[e0], edges[e1]);
            }
        }
    }

    if (segments.empty()) {
        std::cerr << "Validation failed: no contour segments found.\n";
        return 1;
    }

    const double key_scale = 1e6;
    std::unordered_map<Key, std::vector<Key>, KeyHash> adj;
    adj.reserve(segments.size() * 2);

    for (const auto &seg : segments) {
        Key k1 = quantize(seg.first, key_scale);
        Key k2 = quantize(seg.second, key_scale);
        adj[k1].push_back(k2);
        adj[k2].push_back(k1);
    }

    for (const auto &kv : adj) {
        if (kv.second.size() != 2) {
            std::cerr << "Validation failed: contour node degree != 2.\n";
            return 1;
        }
    }

    std::unordered_set<Key, KeyHash> visited;
    std::vector<Key> contour_keys;
    contour_keys.reserve(adj.size());

    double best_perimeter = -1.0;
    for (const auto &kv : adj) {
        const Key start_key = kv.first;
        if (visited.count(start_key)) continue;

        std::vector<Key> loop;
        loop.reserve(1024);
        Key prev{std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::min()};
        Key cur = start_key;
        while (true) {
            if (visited.count(cur)) {
                if (cur == start_key) break;
                std::cerr << "Validation failed: contour loop revisited unexpectedly.\n";
                return 1;
            }
            visited.insert(cur);
            loop.push_back(cur);
            const auto &nbrs = adj[cur];
            Key next = nbrs[0];
            if (prev == next) next = nbrs[1];
            prev = cur;
            cur = next;
            if (loop.size() > adj.size() + 1) {
                std::cerr << "Validation failed: contour walk overflow.\n";
                return 1;
            }
        }

        double per = 0.0;
        for (size_t i = 0; i < loop.size(); ++i) {
            Point a = from_key(loop[i], key_scale);
            Point b = from_key(loop[(i + 1) % loop.size()], key_scale);
            per += dist(a, b);
        }
        if (per > best_perimeter) {
            best_perimeter = per;
            contour_keys = std::move(loop);
        }
    }

    if (contour_keys.empty()) {
        std::cerr << "Validation failed: contour loop not found.\n";
        return 1;
    }

    if (validate && visited.size() != adj.size()) {
        std::cerr << "Validation warning: multiple contour loops detected; using largest.\n";
    }

    std::vector<Point> contour;
    contour.reserve(contour_keys.size());
    for (const auto &k : contour_keys) {
        contour.push_back(from_key(k, key_scale));
    }

    const int m = static_cast<int>(contour.size());
    std::vector<double> prefix(static_cast<size_t>(m) + 1, 0.0);
    for (int i = 0; i < m; ++i) {
        prefix[i + 1] = prefix[i] + dist(contour[i], contour[(i + 1) % m]);
    }
    const double perimeter = prefix[m];

    std::vector<int> visA;
    std::vector<int> visB;
    std::vector<double> distA(m, 0.0);
    std::vector<double> distB(m, 0.0);
    visA.reserve(m);
    visB.reserve(m);

    for (int i = 0; i < m; ++i) {
        distA[i] = dist(A, contour[i]);
        distB[i] = dist(B, contour[i]);
        if (visible(A, contour[i], f_min, vis_step)) visA.push_back(i);
        if (visible(B, contour[i], f_min, vis_step)) visB.push_back(i);
    }

    if (visA.empty() || visB.empty()) {
        std::cerr << "Validation failed: no visible contour points from A or B.\n";
        return 1;
    }

    double best = std::numeric_limits<double>::infinity();
    if (visible(A, B, f_min, vis_step)) best = dist(A, B);

    for (int i : visA) {
        for (int j : visB) {
            double arc = std::abs(prefix[j] - prefix[i]);
            arc = std::min(arc, perimeter - arc);
            double total = distA[i] + distB[j] + arc;
            if (total < best) best = total;
        }
    }

    if (!std::isfinite(best)) {
        std::cerr << "Validation failed: shortest path not found.\n";
        return 1;
    }

    if (validate) {
        const double direct = dist(A, B);
        if (best + 1e-7 < direct) {
            std::cerr << "Validation failed: shortest path shorter than straight line.\n";
            return 1;
        }
        if (contour.size() < 1000) {
            std::cerr << "Validation warning: contour resolution seems too low.\n";
        }
    }

    std::cout << std::fixed << std::setprecision(3) << best << "\n";
    return 0;
}
