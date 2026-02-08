#include <array>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::cout;
using std::int64_t;
using std::size_t;
using std::string;
using std::uint64_t;
using std::vector;

struct Point {
    int x;
    int y;
};

struct DispInfo {
    uint8_t count;
    bool tangent;
    Point p1;
    Point p2;
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

static vector<Point> non_opposite_representatives(const vector<Point>& points) {
    vector<Point> sorted = points;
    std::sort(sorted.begin(), sorted.end(), [](const Point& a, const Point& b) {
        return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    });

    std::unordered_set<int64_t> used;
    used.reserve(points.size() * 2 + 16);

    vector<Point> reps;
    reps.reserve(points.size() / 2 + 2);

    for (const Point& v : sorted) {
        const int64_t k = pair_key(v.x, v.y);
        if (used.find(k) != used.end()) {
            continue;
        }
        const Point ov{-v.x, -v.y};
        used.insert(k);
        used.insert(pair_key(ov.x, ov.y));

        const bool lex_greater = (v.x > ov.x) || (v.x == ov.x && v.y > ov.y);
        reps.push_back(lex_greater ? v : ov);
    }

    return reps;
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

static std::unordered_map<int64_t, DispInfo> precompute_displacements(int m, const vector<Point>& points) {
    std::unordered_set<int64_t> point_set;
    point_set.reserve(points.size() * 2 + 16);

    int max_abs = 0;
    for (const Point& p : points) {
        point_set.insert(pair_key(p.x, p.y));
        max_abs = std::max(max_abs, std::max(std::abs(p.x), std::abs(p.y)));
    }

    std::unordered_map<int64_t, DispInfo> info;
    info.reserve((4 * max_abs + 1) * (4 * max_abs + 1));

    for (int dx = -2 * max_abs; dx <= 2 * max_abs; ++dx) {
        for (int dy = -2 * max_abs; dy <= 2 * max_abs; ++dy) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int64_t d2 = 1LL * dx * dx + 1LL * dy * dy;
            const bool tangent = (d2 == 4LL * m);
            if (d2 > 4LL * m && !tangent) {
                continue;
            }

            int count = 0;
            Point p1{0, 0};
            Point p2{0, 0};

            for (const Point& p : points) {
                if (point_set.find(pair_key(p.x - dx, p.y - dy)) != point_set.end()) {
                    ++count;
                    if (count == 1) {
                        p1 = p;
                    } else if (count == 2) {
                        p2 = p;
                    }
                }
            }

            info.emplace(pair_key(dx, dy), DispInfo{static_cast<uint8_t>(count), tangent, p1, p2});
        }
    }

    return info;
}

static bool test_selected_vectors(
    const vector<Point>& v,
    const std::unordered_map<int64_t, DispInfo>& disp_info,
    const vector<int>& even_masks
) {
    const int d = static_cast<int>(v.size());
    const int full = 1 << d;
    const int circles = 1 << (d - 1);

    vector<Point> sums(full, Point{0, 0});
    for (int mask = 1; mask < full; ++mask) {
        const int b = __builtin_ctz(static_cast<unsigned>(mask));
        const int pm = mask ^ (1 << b);
        sums[mask] = {sums[pm].x + v[b].x, sums[pm].y + v[b].y};
    }

    vector<Point> centers(circles);
    centers.reserve(circles);

    std::unordered_set<int64_t> center_keys;
    center_keys.reserve(circles * 2 + 16);

    for (int i = 0; i < circles; ++i) {
        centers[i] = sums[even_masks[i]];
        const int64_t key = pair_key(centers[i].x, centers[i].y);
        if (center_keys.find(key) != center_keys.end()) {
            return false;
        }
        center_keys.insert(key);
    }

    vector<vector<int>> adj(circles);
    std::unordered_set<int64_t> harmony_keys;
    harmony_keys.reserve(circles * 4 + 16);

    for (int i = 0; i < circles; ++i) {
        const Point& ci = centers[i];
        for (int j = i + 1; j < circles; ++j) {
            const Point& cj = centers[j];
            const int dx = cj.x - ci.x;
            const int dy = cj.y - ci.y;

            const auto it = disp_info.find(pair_key(dx, dy));
            if (it == disp_info.end()) {
                continue;
            }

            const DispInfo& in = it->second;
            if (in.tangent) {
                return false;
            }

            if (in.count == 2) {
                harmony_keys.insert(pair_key(ci.x + in.p1.x, ci.y + in.p1.y));
                harmony_keys.insert(pair_key(ci.x + in.p2.x, ci.y + in.p2.y));

                if (static_cast<int>(harmony_keys.size()) > circles) {
                    return false;
                }

                adj[i].push_back(j);
                adj[j].push_back(i);
            } else if (in.count != 0) {
                return false;
            }
        }
    }

    if (static_cast<int>(harmony_keys.size()) != circles) {
        return false;
    }

    vector<char> visited(circles, 0);
    vector<int> stack;
    stack.reserve(circles);
    stack.push_back(0);
    visited[0] = 1;

    int seen = 0;
    while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        ++seen;

        for (const int vtx : adj[u]) {
            if (!visited[vtx]) {
                visited[vtx] = 1;
                stack.push_back(vtx);
            }
        }
    }

    return seen == circles;
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

    const auto disp_info = precompute_displacements(m, lattice);

    vector<int> even_masks;
    even_masks.reserve(1 << (d - 1));
    for (int mask = 0; mask < (1 << d); ++mask) {
        if ((__builtin_popcount(static_cast<unsigned>(mask)) & 1) == 0) {
            even_masks.push_back(mask);
        }
    }

    const int u = static_cast<int>(pairs.size());
    const int first_max = u - d;
    if (first_max < 0) {
        return false;
    }

    const bool exhaustive_signs = (d <= 6);

    if (thread_count <= 0) {
        thread_count = 1;
    }
    thread_count = std::min(thread_count, first_max + 1);

    std::atomic<bool> found{false};

    auto worker = [&](int tid) {
        vector<int> combo(d, 0);
        vector<Point> selected(d);

        std::function<void(int, int)> dfs = [&](int pos, int next_idx) {
            if (found.load(std::memory_order_relaxed)) {
                return;
            }

            if (pos == d) {
                if (exhaustive_signs) {
                    const int sign_masks = 1 << d;
                    for (int sm = 0; sm < sign_masks; ++sm) {
                        if (found.load(std::memory_order_relaxed)) {
                            return;
                        }
                        if (sm & 1) {
                            continue;
                        }
                        for (int i = 0; i < d; ++i) {
                            selected[i] = pairs[combo[i]][(sm >> i) & 1];
                        }
                        if (test_selected_vectors(selected, disp_info, even_masks)) {
                            found.store(true, std::memory_order_relaxed);
                            return;
                        }
                    }
                } else {
                    for (int i = 0; i < d; ++i) {
                        selected[i] = pairs[combo[i]][0];
                    }
                    if (test_selected_vectors(selected, disp_info, even_masks)) {
                        found.store(true, std::memory_order_relaxed);
                    }
                }
                return;
            }

            const int remaining = d - pos;
            const int limit = u - remaining;
            for (int i = next_idx; i <= limit; ++i) {
                if (found.load(std::memory_order_relaxed)) {
                    return;
                }
                combo[pos] = i;
                dfs(pos + 1, i + 1);
            }
        };

        for (int first = tid; first <= first_max; first += thread_count) {
            if (found.load(std::memory_order_relaxed)) {
                return;
            }
            combo[0] = first;
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

int count_split_primes_distinct(uint64_t m) {
    int count = 0;
    for (uint64_t p = 2; p * p <= m; ++p) {
        if (m % p != 0) {
            continue;
        }
        int exponent = 0;
        while (m % p == 0) {
            m /= p;
            ++exponent;
        }
        if (exponent > 0 && (p % 4 == 1)) {
            ++count;
        }
    }
    if (m > 1 && (m % 4 == 1)) {
        ++count;
    }
    return count;
}

uint64_t product_of_first_split_primes(size_t k) {
    uint64_t product = 1;

    auto is_prime = [](uint64_t x) {
        if (x < 2) {
            return false;
        }
        if (x % 2 == 0) {
            return x == 2;
        }
        for (uint64_t d = 3; d * d <= x; d += 2) {
            if (x % d == 0) {
                return false;
            }
        }
        return true;
    };

    uint64_t p = 2;
    size_t found = 0;
    while (found < k) {
        ++p;
        if (!is_prime(p) || (p % 4 != 1)) {
            continue;
        }

        const __uint128_t tmp = static_cast<__uint128_t>(product) * p;
        if (tmp > std::numeric_limits<uint64_t>::max()) {
            throw std::overflow_error("overflow in product_of_first_split_primes");
        }

        product = static_cast<uint64_t>(tmp);
        ++found;
    }

    return product;
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
        if (find_for_m_and_dimension(m, dimension, thread_count)) {
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
    for (uint64_t n = 2; n <= 8; ++n) {
        const uint64_t cur = solve_r_sq(n, thread_count);
        if (n > 2 && cur < prev) {
            cerr << "Checkpoint failed: monotonicity violated at n=" << n << "\n";
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
