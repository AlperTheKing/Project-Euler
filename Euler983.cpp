#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::cout;
using std::int64_t;
using std::string;
using std::uint64_t;
using std::vector;

struct Point {
    int x;
    int y;
};

static inline int64_t pair_key(int x, int y) {
    return (static_cast<int64_t>(x) << 32) ^ static_cast<uint32_t>(y);
}

static vector<Point> circle_points(int m) {
    const int r = static_cast<int>(std::sqrt(static_cast<double>(m)));
    vector<Point> pts;
    pts.reserve(64);
    for (int x = -r; x <= r; ++x) {
        const int y2 = m - x * x;
        if (y2 < 0) {
            continue;
        }
        const int y = static_cast<int>(std::sqrt(static_cast<double>(y2)));
        if (y * y == y2) {
            pts.push_back({x, y});
            if (y != 0) {
                pts.push_back({x, -y});
            }
        }
    }
    std::sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    });
    pts.erase(std::unique(pts.begin(), pts.end(), [](const Point& a, const Point& b) {
        return a.x == b.x && a.y == b.y;
    }), pts.end());
    return pts;
}

static vector<std::array<Point, 2>> opposite_pairs(const vector<Point>& points) {
    vector<Point> sorted = points;
    std::sort(sorted.begin(), sorted.end(), [](const Point& a, const Point& b) {
        return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    });

    std::unordered_set<int64_t> used;
    used.reserve(points.size() * 2 + 16);

    vector<std::array<Point, 2>> pairs;
    pairs.reserve(points.size() / 2 + 2);

    for (const Point& v : sorted) {
        const int64_t k = pair_key(v.x, v.y);
        if (used.find(k) != used.end()) {
            continue;
        }
        const Point ov{-v.x, -v.y};
        used.insert(k);
        used.insert(pair_key(ov.x, ov.y));
        pairs.push_back({v, ov});
    }

    return pairs;
}

static int opposite_pair_count_by_factorization(uint64_t m) {
    if (m == 0) {
        return 0;
    }

    uint64_t n = m;
    int pairs = 2;

    int e2 = 0;
    while ((n & 1ULL) == 0ULL) {
        n >>= 1;
        ++e2;
    }
    (void)e2;

    for (uint64_t p = 3; p * p <= n; p += 2) {
        if (n % p != 0) {
            continue;
        }
        int exp = 0;
        while (n % p == 0) {
            n /= p;
            ++exp;
        }
        if ((p & 3ULL) == 3ULL) {
            if ((exp & 1) != 0) {
                return 0;
            }
        } else if ((p & 3ULL) == 1ULL) {
            if (pairs > std::numeric_limits<int>::max() / (exp + 1)) {
                return std::numeric_limits<int>::max();
            }
            pairs *= (exp + 1);
        }
    }

    if (n > 1) {
        if ((n & 3ULL) == 3ULL) {
            return 0;
        }
        if ((n & 3ULL) == 1ULL) {
            if (pairs > std::numeric_limits<int>::max() / 2) {
                return std::numeric_limits<int>::max();
            }
            pairs *= 2;
        }
    }

    return pairs;
}

static std::unordered_set<int64_t> build_bad_displacements(const vector<Point>& points) {
    std::unordered_set<int64_t> bad;
    bad.reserve(points.size() * points.size() * 2 + 16);
    for (const Point& a : points) {
        for (const Point& b : points) {
            if (a.x == b.x && a.y == b.y) {
                continue;
            }
            bad.insert(pair_key(a.x - b.x, a.y - b.y));
        }
    }
    return bad;
}

static bool four_tuple_has_forbidden_sum(
    const vector<Point>& selected,
    int used_count,
    const std::unordered_set<int64_t>& bad_disp
) {
    if (used_count < 4) {
        return false;
    }

    const Point& d = selected[used_count - 1];
    for (int i = 0; i < used_count - 1; ++i) {
        const Point& a = selected[i];
        for (int j = i + 1; j < used_count - 1; ++j) {
            const Point& b = selected[j];
            for (int k = j + 1; k < used_count - 1; ++k) {
                const Point& c = selected[k];
                for (int mask = 0; mask < 16; ++mask) {
                    const int sx =
                        ((mask & 1) ? -a.x : a.x) +
                        ((mask & 2) ? -b.x : b.x) +
                        ((mask & 4) ? -c.x : c.x) +
                        ((mask & 8) ? -d.x : d.x);
                    const int sy =
                        ((mask & 1) ? -a.y : a.y) +
                        ((mask & 2) ? -b.y : b.y) +
                        ((mask & 4) ? -c.y : c.y) +
                        ((mask & 8) ? -d.y : d.y);
                    if (sx == 0 && sy == 0) {
                        continue;
                    }
                    if (bad_disp.find(pair_key(sx, sy)) != bad_disp.end()) {
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

static bool test_selected_vectors(
    const vector<Point>& v,
    const vector<Point>& lattice
) {
    const int d = static_cast<int>(v.size());
    const int full = 1 << d;

    vector<Point> sums(full, Point{0, 0});
    for (int mask = 1; mask < full; ++mask) {
        const int bit = __builtin_ctz(static_cast<unsigned>(mask));
        const int pm = mask ^ (1 << bit);
        sums[mask] = {sums[pm].x + v[bit].x, sums[pm].y + v[bit].y};
    }

    std::unordered_set<int64_t> even_keys;
    std::unordered_set<int64_t> odd_keys;
    even_keys.reserve(full + 16);
    odd_keys.reserve(full + 16);
    vector<Point> even_points;
    even_points.reserve(full / 2 + 1);

    for (int mask = 0; mask < full; ++mask) {
        const int64_t key = pair_key(sums[mask].x, sums[mask].y);
        if ((__builtin_popcount(static_cast<unsigned>(mask)) & 1) == 0) {
            if (!even_keys.insert(key).second) {
                return false;
            }
            even_points.push_back(sums[mask]);
        } else {
            if (!odd_keys.insert(key).second) {
                return false;
            }
        }
    }

    for (const Point& c : even_points) {
        for (const Point& dxy : lattice) {
            const int nx = c.x + dxy.x;
            const int ny = c.y + dxy.y;
            if (odd_keys.find(pair_key(nx, ny)) != odd_keys.end()) {
                continue;
            }
            int count = 0;
            for (const Point& p : lattice) {
                if (even_keys.find(pair_key(nx + p.x, ny + p.y)) != even_keys.end()) {
                    ++count;
                    if (count >= 2) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

static bool find_for_m_and_dimension(int m, int d, int thread_count) {
    const vector<Point> lattice = circle_points(m);
    if (lattice.empty()) {
        return false;
    }

    const auto pairs = opposite_pairs(lattice);
    if (static_cast<int>(pairs.size()) < d) {
        return false;
    }

    vector<Point> reps;
    reps.reserve(pairs.size());
    for (const auto& pr : pairs) {
        const Point& a = pr[0];
        const Point& b = pr[1];
        const bool take_a = (a.x > b.x) || (a.x == b.x && a.y > b.y);
        reps.push_back(take_a ? a : b);
    }

    const auto bad_disp = build_bad_displacements(lattice);

    const int u = static_cast<int>(reps.size());
    const int first_max = u - d;
    if (first_max < 0) {
        return false;
    }

    if (thread_count <= 0) {
        thread_count = 1;
    }
    thread_count = std::min(thread_count, first_max + 1);

    std::atomic<bool> found{false};

    auto worker = [&](int tid) {
        vector<Point> selected(d);

        std::function<void(int, int)> dfs = [&](int pos, int next_idx) {
            if (found.load(std::memory_order_relaxed)) {
                return;
            }
            if (pos == d) {
                if (test_selected_vectors(selected, lattice)) {
                    found.store(true, std::memory_order_relaxed);
                }
                return;
            }

            const int remaining = d - pos;
            const int limit = u - remaining;
            for (int i = next_idx; i <= limit; ++i) {
                if (found.load(std::memory_order_relaxed)) {
                    return;
                }
                selected[pos] = reps[i];
                if (four_tuple_has_forbidden_sum(selected, pos + 1, bad_disp)) {
                    continue;
                }
                dfs(pos + 1, i + 1);
            }
        };

        for (int first = tid; first <= first_max; first += thread_count) {
            if (found.load(std::memory_order_relaxed)) {
                return;
            }
            selected[0] = reps[first];
            dfs(1, first + 1);
        }
    };

    vector<std::thread> threads;
    threads.reserve(thread_count);
    for (int t = 0; t < thread_count; ++t) {
        threads.emplace_back(worker, t);
    }
    for (auto& th : threads) {
        th.join();
    }

    return found.load(std::memory_order_relaxed);
}

uint64_t solve_r_sq(uint64_t n, int thread_count = 1) {
    if (n <= 2) {
        return 1;
    }

    int d = 1;
    uint64_t capacity = 1;
    while (capacity < n) {
        ++d;
        if (d >= 64) {
            capacity = std::numeric_limits<uint64_t>::max();
            break;
        }
        capacity <<= 1;
    }

    const int dimension = d;

    int m = 1;
    while (true) {
        if (opposite_pair_count_by_factorization(static_cast<uint64_t>(m)) >= dimension &&
            find_for_m_and_dimension(m, dimension, thread_count)) {
            return static_cast<uint64_t>(m);
        }
        ++m;
    }
}

bool run_checkpoints(int thread_count) {
    const uint64_t r2 = solve_r_sq(2, thread_count);
    if (r2 != 1) {
        cerr << "Checkpoint failed: expected R(2)^2 = 1, got " << r2 << "\n";
        return false;
    }

    const uint64_t r4 = solve_r_sq(4, thread_count);
    if (r4 != 5) {
        cerr << "Checkpoint failed: expected R(4)^2 = 5, got " << r4 << "\n";
        return false;
    }

    const uint64_t r4_again = solve_r_sq(4, thread_count);
    if (r4_again != r4) {
        cerr << "Checkpoint failed: recomputation mismatch for n=4 ("
             << r4 << " vs " << r4_again << ")\n";
        return false;
    }

    uint64_t prev = 0;
    for (uint64_t ncheck = 2; ncheck <= 8; ++ncheck) {
        const uint64_t cur = solve_r_sq(ncheck, thread_count);
        if (ncheck > 2 && cur < prev) {
            cerr << "Checkpoint failed: monotonicity violated at n=" << ncheck << "\n";
            return false;
        }
        prev = cur;
    }

    return true;
}

int main(int argc, char** argv) {
    uint64_t n = 500;
    bool skip_checkpoints = false;

    int thread_count = static_cast<int>(std::thread::hardware_concurrency());
    if (thread_count <= 0) {
        thread_count = 1;
    }

    for (int i = 1; i < argc; ++i) {
        const string arg(argv[i]);
        if (arg.rfind("--n=", 0) == 0) {
            n = std::stoull(arg.substr(4));
        } else if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else if (arg.rfind("--threads=", 0) == 0) {
            thread_count = std::stoi(arg.substr(10));
            if (thread_count <= 0) {
                thread_count = 1;
            }
        } else {
            cerr << "Unknown argument: " << arg << "\n";
            return 1;
        }
    }

    if (!skip_checkpoints) {
        if (!run_checkpoints(thread_count)) {
            return 1;
        }
    }

    cout << solve_r_sq(n, thread_count) << "\n";
    return 0;
}
