#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <thread>
#include <vector>

namespace {

struct Point {
    int x;
    int y;
    int z;
};

constexpr long double kPi = 3.141592653589793238462643383279502884L;

long double edge_risk(const Point& a, const Point& b, int r) {
    const long double dot = static_cast<long double>(a.x) * b.x +
                            static_cast<long double>(a.y) * b.y +
                            static_cast<long double>(a.z) * b.z;
    long double c = dot / (static_cast<long double>(r) * static_cast<long double>(r));
    if (c > 1.0L) c = 1.0L;
    if (c < -1.0L) c = -1.0L;
    const long double t = std::acos(c) / kPi;
    return t * t;
}

std::vector<Point> generate_all_points(int r) {
    const int rr = r * r;
    std::vector<Point> pts;
    pts.reserve(220000);
    for (int x = -r; x <= r; ++x) {
        const int x2 = x * x;
        for (int y = -r; y <= r; ++y) {
            const int z2 = rr - x2 - y * y;
            if (z2 < 0) continue;
            const int z = static_cast<int>(std::llround(std::sqrt(static_cast<long double>(z2))));
            if (z * z != z2) continue;
            pts.push_back({x, y, z});
            if (z != 0) pts.push_back({x, y, -z});
        }
    }
    return pts;
}

std::pair<int, int> pole_indices(const std::vector<Point>& pts, int r) {
    int north = -1;
    int south = -1;
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
        if (pts[static_cast<std::size_t>(i)].x == 0 && pts[static_cast<std::size_t>(i)].y == 0 &&
            pts[static_cast<std::size_t>(i)].z == r) {
            north = i;
        }
        if (pts[static_cast<std::size_t>(i)].x == 0 && pts[static_cast<std::size_t>(i)].y == 0 &&
            pts[static_cast<std::size_t>(i)].z == -r) {
            south = i;
        }
    }
    return {north, south};
}

long double dijkstra_complete(const std::vector<Point>& pts, int src, int dst, int r) {
    const int n = static_cast<int>(pts.size());
    std::vector<long double> dist(static_cast<std::size_t>(n),
                                  std::numeric_limits<long double>::infinity());
    std::vector<unsigned char> done(static_cast<std::size_t>(n), 0);
    dist[static_cast<std::size_t>(src)] = 0.0L;
    int u = src;
    while (true) {
        done[static_cast<std::size_t>(u)] = 1;
        if (u == dst) break;
        const long double du = dist[static_cast<std::size_t>(u)];
        for (int v = 0; v < n; ++v) {
            if (done[static_cast<std::size_t>(v)] || v == u) continue;
            const long double nd = du + edge_risk(pts[static_cast<std::size_t>(u)],
                                                  pts[static_cast<std::size_t>(v)], r);
            if (nd < dist[static_cast<std::size_t>(v)]) dist[static_cast<std::size_t>(v)] = nd;
        }
        int best = -1;
        long double bd = std::numeric_limits<long double>::infinity();
        for (int v = 0; v < n; ++v) {
            if (done[static_cast<std::size_t>(v)]) continue;
            if (dist[static_cast<std::size_t>(v)] < bd) {
                bd = dist[static_cast<std::size_t>(v)];
                best = v;
            }
        }
        if (best < 0) break;
        u = best;
    }
    return dist[static_cast<std::size_t>(dst)];
}

std::vector<Point> generate_reduced_points(int r, std::vector<long double>& asin_pi) {
    asin_pi.resize(static_cast<std::size_t>(r + 1));
    const long double inv_r = 1.0L / static_cast<long double>(r);
    for (int z = 0; z <= r; ++z) {
        asin_pi[static_cast<std::size_t>(z)] = std::asin(static_cast<long double>(z) * inv_r) / kPi;
    }

    std::vector<Point> pts;
    pts.reserve(30000);
    const int rr = r * r;

    int z0 = r;
    for (int x = 0; 3LL * x * x < rr; ++x) {
        int y = x;
        int z = z0;
        long long h = 1LL * x * x + 1LL * y * y + 1LL * z * z - rr;
        while (h > 0) {
            h -= 2LL * z - 1LL;
            --z;
        }
        z0 = z;

        while (y <= z) {
            if (h == 0) {
                if (y == 0) {
                    pts.push_back({0, 0, r});
                    pts.push_back({0, r, 0});
                } else if (x == 0) {
                    pts.push_back({0, y, z});
                    pts.push_back({0, z, y});
                    pts.push_back({y, z, 0});
                } else if (y == z) {
                    pts.push_back({x, y, y});
                    pts.push_back({y, y, x});
                } else if (x == y) {
                    pts.push_back({x, x, z});
                    pts.push_back({x, z, x});
                } else {
                    pts.push_back({x, y, z});
                    pts.push_back({x, z, y});
                    pts.push_back({y, z, x});
                }
            }

            h += 2LL * y + 1LL;
            ++y;
            if (h > 0) {
                h -= 2LL * z - 1LL;
                --z;
            }
        }
    }

    std::sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        if (a.z != b.z) return a.z > b.z;
        if (a.y != b.y) return a.y > b.y;
        return a.x > b.x;
    });
    pts.erase(std::unique(pts.begin(), pts.end(),
                          [](const Point& a, const Point& b) {
                              return a.x == b.x && a.y == b.y && a.z == b.z;
                          }),
              pts.end());
    return pts;
}

long double dijkstra_pruned(const std::vector<Point>& pts,
                            const std::vector<long double>& asin_pi,
                            int r,
                            long double riskmax) {
    const int n = static_cast<int>(pts.size());
    if (n < 2) return std::numeric_limits<long double>::infinity();

    std::vector<long double> risk(static_cast<std::size_t>(n), 2.0L);
    std::vector<unsigned char> done(static_cast<std::size_t>(n), 0);

    int cur = 0;
    risk[0] = 0.0L;
    while (true) {
        done[static_cast<std::size_t>(cur)] = 1;
        const Point& pc = pts[static_cast<std::size_t>(cur)];
        const long double rc = risk[static_cast<std::size_t>(cur)];

        if (riskmax > 0.0L) {
            const long double dmax = riskmax - 2.0L * rc;
            for (int j = 0; j < n; ++j) {
                if (done[static_cast<std::size_t>(j)]) continue;
                const Point& pj = pts[static_cast<std::size_t>(j)];
                long double f = asin_pi[static_cast<std::size_t>(pc.z)] -
                                asin_pi[static_cast<std::size_t>(pj.z)];
                if (f * f > dmax) {
                    if (j > cur) break;
                    continue;
                }
                f = asin_pi[static_cast<std::size_t>(pc.x)] -
                    asin_pi[static_cast<std::size_t>(pj.x)];
                if (f * f > dmax) continue;
                f = asin_pi[static_cast<std::size_t>(pc.y)] -
                    asin_pi[static_cast<std::size_t>(pj.y)];
                if (f * f > dmax) continue;

                const long double nd = rc + edge_risk(pc, pj, r);
                if (nd < risk[static_cast<std::size_t>(j)]) {
                    risk[static_cast<std::size_t>(j)] = nd;
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                if (done[static_cast<std::size_t>(j)] || j == cur) continue;
                const long double nd = rc + edge_risk(pc, pts[static_cast<std::size_t>(j)], r);
                if (nd < risk[static_cast<std::size_t>(j)]) {
                    risk[static_cast<std::size_t>(j)] = nd;
                }
            }
        }

        int next = -1;
        long double best = 2.0L;
        for (int j = 0; j < n; ++j) {
            if (done[static_cast<std::size_t>(j)]) continue;
            if (risk[static_cast<std::size_t>(j)] < best) {
                best = risk[static_cast<std::size_t>(j)];
                next = j;
            }
        }
        if (next < 0) break;
        cur = next;
    }

    long double ans = 2.0L;
    for (int i = 0; i < n; ++i) {
        const long double e = 2.0L * asin_pi[static_cast<std::size_t>(pts[static_cast<std::size_t>(i)].z)];
        const long double cand = 2.0L * risk[static_cast<std::size_t>(i)] + e * e;
        if (cand < ans) ans = cand;
    }
    return ans;
}

std::vector<Point> downsample(const std::vector<Point>& pts, int stride) {
    std::vector<Point> out;
    out.reserve(pts.size() / static_cast<std::size_t>(stride) + 2);
    for (std::size_t i = 0; i < pts.size(); i += static_cast<std::size_t>(stride)) {
        out.push_back(pts[i]);
    }
    return out;
}

long double minimal_risk_fast(int r) {
    std::vector<long double> asin_pi;
    const std::vector<Point> pts = generate_reduced_points(r, asin_pi);
    if (pts.size() < 2) return std::numeric_limits<long double>::infinity();

    long double d = 0.0L;
    for (int stride : {27, 9, 3, 1}) {
        std::vector<Point> sample = (stride == 1) ? pts : downsample(pts, stride);
        if (sample.size() < 2) continue;
        d = dijkstra_pruned(sample, asin_pi, r, d);
    }
    return d;
}

bool approx_equal_10dp(long double a, long double b) {
    return std::fabsl(a - b) < 0.5e-10L;
}

bool run_checkpoints() {
    const long double m7 = minimal_risk_fast(7);
    if (!approx_equal_10dp(m7, 0.1784943998L)) {
        std::cerr << "Checkpoint failed: M(7)\n";
        return false;
    }

    const int r = 31;
    const auto pts = generate_all_points(r);
    const auto [north, south] = pole_indices(pts, r);
    const long double exact = dijkstra_complete(pts, north, south, r);
    const long double fast = minimal_risk_fast(r);
    if (std::fabsl(exact - fast) > 1e-12L) {
        std::cerr << "Checkpoint failed: M(31) exact mismatch\n";
        return false;
    }

    return true;
}

long double solve() {
    struct Task {
        std::atomic<int>* next_n;
        std::array<long double, 16>* values;
    };

    auto worker = [](void* arg) -> void* {
        auto* task = static_cast<Task*>(arg);
        while (true) {
            const int n = task->next_n->fetch_add(1);
            if (n > 15) break;
            const int r = (1 << n) - 1;
            (*task->values)[static_cast<std::size_t>(n)] = minimal_risk_fast(r);
        }
        return nullptr;
    };

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) hw = 1U;
    const unsigned thread_count = std::min<unsigned>(8U, std::min<unsigned>(15U, hw));

    std::array<long double, 16> values{};
    std::atomic<int> next_n{1};
    Task task{&next_n, &values};
    std::vector<pthread_t> threads(thread_count);
    for (unsigned t = 0; t < thread_count; ++t) {
        pthread_create(&threads[static_cast<std::size_t>(t)], nullptr, worker, &task);
    }
    for (unsigned t = 0; t < thread_count; ++t) {
        pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
    }

    long double sum = 0.0L;
    for (int n = 1; n <= 15; ++n) {
        sum += values[static_cast<std::size_t>(n)];
    }
    return sum;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
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
