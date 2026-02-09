#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Point {
    int x;
    int y;
    int z;
};

int floor_div(const int a, const int b) {
    if (a >= 0) {
        return a / b;
    }
    return -(((-a) + b - 1) / b);
}

std::uint64_t pack_cell_key(const int cx, const int cy, const int cz) {
    // Cell coordinates stay small for our ranges; shift into nonnegative space.
    const std::uint64_t ox = static_cast<std::uint64_t>(cx + 2048);
    const std::uint64_t oy = static_cast<std::uint64_t>(cy + 2048);
    const std::uint64_t oz = static_cast<std::uint64_t>(cz + 2048);
    return (ox << 24U) ^ (oy << 12U) ^ oz;
}

std::vector<Point> generate_points(const int r) {
    const int rr = r * r;
    std::vector<Point> pts;
    pts.reserve(200000);

    for (int x = -r; x <= r; ++x) {
        const int x2 = x * x;
        for (int y = -r; y <= r; ++y) {
            const int z2 = rr - x2 - y * y;
            if (z2 < 0) {
                continue;
            }
            const int z = static_cast<int>(std::llround(std::sqrt(static_cast<long double>(z2))));
            if (z * z != z2) {
                continue;
            }
            pts.push_back({x, y, z});
            if (z != 0) {
                pts.push_back({x, y, -z});
            }
        }
    }

    return pts;
}

long double edge_risk(const Point& a, const Point& b, const long double inv_rr) {
    const long double dot = static_cast<long double>(a.x) * b.x +
                            static_cast<long double>(a.y) * b.y +
                            static_cast<long double>(a.z) * b.z;
    long double c = dot * inv_rr;
    if (c > 1.0L) {
        c = 1.0L;
    } else if (c < -1.0L) {
        c = -1.0L;
    }
    const long double theta = std::acos(c);
    const long double x = theta / std::acos(-1.0L);
    return x * x;
}

long double dijkstra_complete(const std::vector<Point>& pts, const int src, const int dst, const int r) {
    const int n = static_cast<int>(pts.size());
    const long double inv_rr = 1.0L / static_cast<long double>(r) / static_cast<long double>(r);

    std::vector<long double> dist(static_cast<std::size_t>(n), std::numeric_limits<long double>::infinity());
    using Node = std::pair<long double, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    dist[static_cast<std::size_t>(src)] = 0.0L;
    pq.push({0.0L, src});

    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[static_cast<std::size_t>(u)]) {
            continue;
        }
        if (u == dst) {
            return d;
        }

        for (int v = 0; v < n; ++v) {
            if (v == u) {
                continue;
            }
            const long double nd = d + edge_risk(pts[static_cast<std::size_t>(u)], pts[static_cast<std::size_t>(v)], inv_rr);
            if (nd < dist[static_cast<std::size_t>(v)]) {
                dist[static_cast<std::size_t>(v)] = nd;
                pq.push({nd, v});
            }
        }
    }

    return dist[static_cast<std::size_t>(dst)];
}

long double dijkstra_sparse(const std::vector<Point>& pts, const int src, const int dst, const int r, const int chord_limit) {
    const int n = static_cast<int>(pts.size());
    const int D = std::max(1, chord_limit);
    const int D2 = D * D;
    const long double inv_rr = 1.0L / static_cast<long double>(r) / static_cast<long double>(r);

    std::unordered_map<std::uint64_t, std::vector<int>> bins;
    bins.reserve(static_cast<std::size_t>(n * 2));

    std::vector<std::array<int, 3>> cell_of(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const Point& p = pts[static_cast<std::size_t>(i)];
        const int cx = floor_div(p.x, D);
        const int cy = floor_div(p.y, D);
        const int cz = floor_div(p.z, D);
        cell_of[static_cast<std::size_t>(i)] = {cx, cy, cz};
        bins[pack_cell_key(cx, cy, cz)].push_back(i);
    }

    std::vector<std::vector<std::pair<int, long double>>> adj(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) {
        const auto c = cell_of[static_cast<std::size_t>(i)];
        const Point& p = pts[static_cast<std::size_t>(i)];

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto it = bins.find(pack_cell_key(c[0] + dx, c[1] + dy, c[2] + dz));
                    if (it == bins.end()) {
                        continue;
                    }
                    for (int j : it->second) {
                        if (j <= i) {
                            continue;
                        }
                        const Point& q = pts[static_cast<std::size_t>(j)];
                        const int ddx = p.x - q.x;
                        const int ddy = p.y - q.y;
                        const int ddz = p.z - q.z;
                        const int chord2 = ddx * ddx + ddy * ddy + ddz * ddz;
                        if (chord2 > D2) {
                            continue;
                        }

                        const long double w = edge_risk(p, q, inv_rr);
                        adj[static_cast<std::size_t>(i)].push_back({j, w});
                        adj[static_cast<std::size_t>(j)].push_back({i, w});
                    }
                }
            }
        }
    }

    std::vector<long double> dist(static_cast<std::size_t>(n), std::numeric_limits<long double>::infinity());
    using Node = std::pair<long double, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    dist[static_cast<std::size_t>(src)] = 0.0L;
    pq.push({0.0L, src});

    while (!pq.empty()) {
        const auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[static_cast<std::size_t>(u)]) {
            continue;
        }
        if (u == dst) {
            return d;
        }

        for (const auto& [v, w] : adj[static_cast<std::size_t>(u)]) {
            const long double nd = d + w;
            if (nd < dist[static_cast<std::size_t>(v)]) {
                dist[static_cast<std::size_t>(v)] = nd;
                pq.push({nd, v});
            }
        }
    }

    return std::numeric_limits<long double>::infinity();
}

std::pair<int, int> pole_indices(const std::vector<Point>& pts, const int r) {
    int north = -1;
    int south = -1;
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
        const Point& p = pts[static_cast<std::size_t>(i)];
        if (p.x == 0 && p.y == 0 && p.z == r) {
            north = i;
        }
        if (p.x == 0 && p.y == 0 && p.z == -r) {
            south = i;
        }
    }
    return {north, south};
}

long double minimal_risk(const int r, const bool force_complete_small = true) {
    const std::vector<Point> pts = generate_points(r);
    const auto [north, south] = pole_indices(pts, r);

    if (north < 0 || south < 0) {
        return std::numeric_limits<long double>::infinity();
    }

    if (force_complete_small && static_cast<int>(pts.size()) <= 4500) {
        return dijkstra_complete(pts, north, south, r);
    }

    const long double root_r = std::sqrt(static_cast<long double>(r));
    int D = static_cast<int>(std::ceil(3.1L * root_r));

    for (int attempt = 0; attempt < 8; ++attempt) {
        const long double risk = dijkstra_sparse(pts, north, south, r, D);
        if (std::isfinite(risk)) {
            return risk;
        }
        D = static_cast<int>(std::ceil(D * 1.25L));
    }

    return std::numeric_limits<long double>::infinity();
}

bool approx_equal_10dp(const long double a, const long double b) {
    return std::fabsl(a - b) < 0.5e-10L;
}

bool run_checkpoints() {
    const long double m7 = minimal_risk(7, true);
    if (!approx_equal_10dp(m7, 0.1784943998L)) {
        std::cerr << "Checkpoint failed: M(7)\n";
        return false;
    }

    // Additional internal consistency checks against exact complete-graph runs.
    for (int r : {31, 63, 127}) {
        const std::vector<Point> pts = generate_points(r);
        const auto [north, south] = pole_indices(pts, r);
        const long double exact = dijkstra_complete(pts, north, south, r);
        const int D = static_cast<int>(std::ceil(3.1L * std::sqrt(static_cast<long double>(r))));
        const long double sparse = dijkstra_sparse(pts, north, south, r, D);
        if (!std::isfinite(sparse) || std::fabsl(exact - sparse) > 1e-12L) {
            std::cerr << "Checkpoint failed: sparse/complete mismatch for r=" << r << '\n';
            return false;
        }
    }

    return true;
}

long double solve() {
    long double sum = 0.0L;
    for (int n = 1; n <= 15; ++n) {
        const int r = (1 << n) - 1;
        sum += minimal_risk(r, false);
    }
    return sum;
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

    const long double answer = solve();
    std::cout << std::fixed << std::setprecision(10) << answer << '\n';
    return 0;
}
