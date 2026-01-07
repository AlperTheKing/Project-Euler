#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;

using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

struct Point {
    i64 x;
    i64 y;
};

struct Dir {
    i64 x;
    i64 y;
    bool operator==(const Dir& other) const { return x == other.x && y == other.y; }
};

struct DirHash {
    size_t operator()(const Dir& d) const {
        u64 hx = static_cast<u64>(d.x) ^ (static_cast<u64>(d.x) >> 32);
        u64 hy = static_cast<u64>(d.y) ^ (static_cast<u64>(d.y) >> 32);
        return static_cast<size_t>((hx * 0x9e3779b97f4a7c15ULL) ^ hy);
    }
};

struct PairEntry {
    uint32_t ab; // packed endpoints
    u64 len;
};

struct Group {
    Dir key;
    vector<PairEntry> pairs;
};

struct Divisor {
    u64 d;
    array<uint8_t, 8> exps;
    uint32_t n_points;
};

static inline i64 llabs64(i64 v) {
    return v < 0 ? -v : v;
}

Point multiply(const Point& a, const Point& b) {
    __int128 x = static_cast<__int128>(a.x) * b.x - static_cast<__int128>(a.y) * b.y;
    __int128 y = static_cast<__int128>(a.x) * b.y + static_cast<__int128>(a.y) * b.x;
    return {static_cast<i64>(x), static_cast<i64>(y)};
}

Point canonical(const Point& z) {
    Point cands[4] = {
        {z.x, z.y},
        {-z.y, z.x},
        {-z.x, -z.y},
        {z.y, -z.x}
    };
    for (const auto& c : cands) {
        if (c.x > 0 && c.y >= 0) return c;
        if (c.x == 0 && c.y > 0) return c;
    }
    return cands[0];
}

bool point_less(const Point& a, const Point& b) {
    if (a.x != b.x) return a.x < b.x;
    return a.y < b.y;
}

Dir normalize_dir(i64 x, i64 y) {
    i64 ax = llabs64(x);
    i64 ay = llabs64(y);
    i64 g = std::gcd(ax, ay);
    x /= g;
    y /= g;
    if (x < 0 || (x == 0 && y < 0)) {
        x = -x;
        y = -y;
    }
    return {x, y};
}

Dir perpendicular_dir(const Dir& d) {
    return normalize_dir(-d.y, d.x);
}

Point find_gaussian_prime(int p) {
    for (int a = 1; a * a <= p; ++a) {
        int b2 = p - a * a;
        int b = static_cast<int>(std::sqrt(b2));
        if (b * b == b2) {
            return {a, b};
        }
    }
    return {0, 0};
}

vector<vector<vector<Point>>> precompute_reps(const vector<int>& primes, const vector<int>& max_exp) {
    vector<vector<vector<Point>>> reps(primes.size());
    for (size_t i = 0; i < primes.size(); ++i) {
        int p = primes[i];
        int emax = max_exp[i];
        Point g = find_gaussian_prime(p);
        Point gc = {g.x, -g.y};

        vector<Point> pow_g(emax + 1), pow_gc(emax + 1);
        pow_g[0] = {1, 0};
        pow_gc[0] = {1, 0};
        for (int e = 1; e <= emax; ++e) {
            pow_g[e] = multiply(pow_g[e - 1], g);
            pow_gc[e] = multiply(pow_gc[e - 1], gc);
        }

        reps[i].resize(emax + 1);
        for (int e = 0; e <= emax; ++e) {
            vector<Point> list;
            list.reserve(e + 1);
            for (int k = 0; k <= e; ++k) {
                Point z = multiply(pow_g[k], pow_gc[e - k]);
                list.push_back(canonical(z));
            }
            sort(list.begin(), list.end(), point_less);
            list.erase(unique(list.begin(), list.end(), [](const Point& a, const Point& b) {
                return a.x == b.x && a.y == b.y;
            }), list.end());
            reps[i][e] = std::move(list);
        }
    }
    return reps;
}

vector<Point> build_points(const array<uint8_t, 8>& exps,
                           const vector<vector<vector<Point>>>& reps) {
    vector<Point> base;
    base.push_back({1, 0});
    for (size_t i = 0; i < reps.size(); ++i) {
        const auto& list = reps[i][exps[i]];
        vector<Point> next;
        next.reserve(base.size() * list.size());
        for (const auto& z : base) {
            for (const auto& r : list) {
                next.push_back(multiply(z, r));
            }
        }
        base.swap(next);
    }

    for (auto& z : base) {
        z = canonical(z);
    }
    sort(base.begin(), base.end(), point_less);
    base.erase(unique(base.begin(), base.end(), [](const Point& a, const Point& b) {
        return a.x == b.x && a.y == b.y;
    }), base.end());

    vector<Point> points;
    points.reserve(base.size() * 4);
    for (const auto& z : base) {
        points.push_back({z.x, z.y});
        points.push_back({-z.y, z.x});
        points.push_back({-z.x, -z.y});
        points.push_back({z.y, -z.x});
    }
    sort(points.begin(), points.end(), point_less);
    points.erase(unique(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x == b.x && a.y == b.y;
    }), points.end());

    return points;
}

struct OverlapScratch {
    vector<uint32_t> degA;
    vector<uint32_t> degB;
    vector<size_t> offA;
    vector<size_t> offB;
    vector<size_t> curA;
    vector<size_t> curB;
    vector<u64> bufA;
    vector<u64> bufB;

    void reset(size_t n) {
        if (degA.size() != n) {
            degA.assign(n, 0);
            degB.assign(n, 0);
        } else {
            fill(degA.begin(), degA.end(), 0);
            fill(degB.begin(), degB.end(), 0);
        }
        if (offA.size() != n + 1) {
            offA.assign(n + 1, 0);
            offB.assign(n + 1, 0);
            curA.resize(n + 1);
            curB.resize(n + 1);
        }
    }
};

u64 count_overlaps(const vector<PairEntry>& A, const vector<PairEntry>& B, size_t n, u64 limit,
                   OverlapScratch& scratch) {
    scratch.reset(n);
    for (const auto& p : A) {
        uint32_t a = p.ab & 0xFFFFu;
        uint32_t b = p.ab >> 16;
        scratch.degA[a]++;
        scratch.degA[b]++;
    }
    for (const auto& p : B) {
        uint32_t a = p.ab & 0xFFFFu;
        uint32_t b = p.ab >> 16;
        scratch.degB[a]++;
        scratch.degB[b]++;
    }

    scratch.offA[0] = 0;
    scratch.offB[0] = 0;
    for (size_t i = 0; i < n; ++i) {
        scratch.offA[i + 1] = scratch.offA[i] + scratch.degA[i];
        scratch.offB[i + 1] = scratch.offB[i] + scratch.degB[i];
    }

    size_t totalA = scratch.offA[n];
    size_t totalB = scratch.offB[n];
    scratch.bufA.resize(totalA);
    scratch.bufB.resize(totalB);
    scratch.curA = scratch.offA;
    scratch.curB = scratch.offB;

    for (const auto& p : A) {
        uint32_t a = p.ab & 0xFFFFu;
        uint32_t b = p.ab >> 16;
        scratch.bufA[scratch.curA[a]++] = p.len;
        scratch.bufA[scratch.curA[b]++] = p.len;
    }
    for (const auto& p : B) {
        uint32_t a = p.ab & 0xFFFFu;
        uint32_t b = p.ab >> 16;
        scratch.bufB[scratch.curB[a]++] = p.len;
        scratch.bufB[scratch.curB[b]++] = p.len;
    }

    for (size_t i = 0; i < n; ++i) {
        size_t a0 = scratch.offA[i];
        size_t a1 = scratch.offA[i + 1];
        size_t b0 = scratch.offB[i];
        size_t b1 = scratch.offB[i + 1];
        if (a1 - a0 > 1) {
            sort(scratch.bufA.begin() + a0, scratch.bufA.begin() + a1);
        }
        if (b1 - b0 > 1) {
            sort(scratch.bufB.begin() + b0, scratch.bufB.begin() + b1);
        }
    }

    // Count pairs that share an endpoint (one point) and satisfy the length limit.
    u64 overlap = 0;
    for (size_t i = 0; i < n; ++i) {
        size_t a0 = scratch.offA[i];
        size_t a1 = scratch.offA[i + 1];
        size_t b0 = scratch.offB[i];
        size_t b1 = scratch.offB[i + 1];
        if (a0 == a1 || b0 == b1) continue;
        size_t j = b1;
        for (size_t k = a0; k < a1; ++k) {
            while (j > b0) {
                u128 sum = static_cast<u128>(scratch.bufA[k]) + scratch.bufB[j - 1];
                if (sum <= limit) break;
                --j;
            }
            if (j == b0) break;
            overlap += static_cast<u64>(j - b0);
        }
    }

    return overlap;
}

u64 count_nonzero(const vector<Point>& pts, u64 d) {
    size_t n = pts.size();
    if (n < 4) return 0;

    unordered_map<Dir, size_t, DirHash> index;
    vector<Group> groups;
    index.reserve(n * 8);
    groups.reserve(n * 8);

    auto get_group = [&](const Dir& key) -> size_t {
        auto it = index.find(key);
        if (it != index.end()) return it->second;
        size_t idx = groups.size();
        index.emplace(key, idx);
        groups.push_back({key, {}});
        return idx;
    };

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            i64 sx = pts[i].x + pts[j].x;
            i64 sy = pts[i].y + pts[j].y;
            if (sx == 0 && sy == 0) continue;
            Dir key = normalize_dir(sx, sy);
            __int128 lx = static_cast<__int128>(sx) * sx;
            __int128 ly = static_cast<__int128>(sy) * sy;
            u64 len = static_cast<u64>(lx + ly);
            size_t idx = get_group(key);
            uint32_t ab = static_cast<uint32_t>(i) | (static_cast<uint32_t>(j) << 16);
            groups[idx].pairs.push_back({ab, len});
        }
    }

    for (auto& g : groups) {
        sort(g.pairs.begin(), g.pairs.end(), [](const PairEntry& a, const PairEntry& b) {
            return a.len < b.len;
        });
    }

    // For nonzero sums, the opposite-pair chords cross iff |u|^2 + |v|^2 <= 4d.
    u64 limit = 4 * d;
    u64 total = 0;
    vector<char> done(groups.size(), 0);
    OverlapScratch scratch;

    for (size_t i = 0; i < groups.size(); ++i) {
        if (done[i]) continue;
        Dir key = groups[i].key;
        Dir perp = perpendicular_dir(key);
        auto it = index.find(perp);
        done[i] = 1;
        if (it == index.end()) continue;
        size_t j = it->second;
        if (done[j]) continue;
        done[j] = 1;

        const auto& A = groups[i].pairs;
        const auto& B = groups[j].pairs;
        if (A.empty() || B.empty()) continue;

        u64 cnt = 0;
        size_t jj = B.size();
        for (const auto& pa : A) {
            while (jj > 0) {
                u128 sum = static_cast<u128>(pa.len) + B[jj - 1].len;
                if (sum <= limit) break;
                --jj;
            }
            if (jj == 0) break;
            cnt += static_cast<u64>(jj);
        }

        u64 overlap = count_overlaps(A, B, n, limit, scratch);
        total += cnt - overlap;
    }

    return total;
}

u64 compute_f(u64 d, const array<uint8_t, 8>& exps,
             const vector<vector<vector<Point>>>& reps) {
    auto points = build_points(exps, reps);
    size_t n = points.size();
    if (n < 4) return 0;

    u64 nonzero = count_nonzero(points, d);

    // One antipodal pair (sum vector 0) is always perpendicular; count those combinatorially.
    u64 m = static_cast<u64>(n / 2);
    u64 zero = m * (m - 1) * (2 * m - 3) / 2;
    return zero + nonzero;
}

void gen_divisors(size_t idx, u64 cur, const vector<int>& primes, const vector<int>& max_exp,
                  array<uint8_t, 8>& cur_exp, vector<Divisor>& out) {
    if (idx == primes.size()) {
        u64 n_points = 4;
        for (size_t i = 0; i < primes.size(); ++i) {
            n_points *= static_cast<u64>(cur_exp[i] + 1);
        }
        out.push_back({cur, cur_exp, static_cast<uint32_t>(n_points)});
        return;
    }
    u64 val = 1;
    for (int e = 0; e <= max_exp[idx]; ++e) {
        cur_exp[idx] = static_cast<uint8_t>(e);
        gen_divisors(idx + 1, cur * val, primes, max_exp, cur_exp, out);
        val *= static_cast<u64>(primes[idx]);
    }
}

u64 compute_S(const vector<int>& primes, const vector<int>& max_exp,
              const vector<vector<vector<Point>>>& reps, bool use_threads, bool show_progress) {
    vector<Divisor> divisors;
    array<uint8_t, 8> cur{};
    gen_divisors(0, 1, primes, max_exp, cur, divisors);

    uint32_t max_n = 0;
    for (const auto& div : divisors) {
        max_n = max(max_n, div.n_points);
    }

    sort(divisors.begin(), divisors.end(), [](const Divisor& a, const Divisor& b) {
        return a.n_points > b.n_points;
    });

    unsigned int threads = 1;
    if (use_threads) {
        unsigned int hw = thread::hardware_concurrency();
        if (hw == 0) hw = 2;
        threads = min(hw, 4u);
        if (max_n > 8000) threads = min(threads, 2u);
    }

    vector<u64> results(divisors.size(), 0);
    atomic<size_t> next(0);
    atomic<size_t> done(0);

    auto worker = [&]() {
        while (true) {
            size_t idx = next.fetch_add(1);
            if (idx >= divisors.size()) break;
            const auto& div = divisors[idx];
            results[idx] = compute_f(div.d, div.exps, reps);
            done.fetch_add(1, std::memory_order_relaxed);
        }
    };

    vector<thread> pool;
    pool.reserve(threads);
    thread monitor;
    if (show_progress) {
        monitor = thread([&]() {
            size_t total = divisors.size();
            while (true) {
                size_t finished = done.load(std::memory_order_relaxed);
                u64 percent = total ? static_cast<u64>(finished * 100 / total) : 100;
                cerr << '\r' << "Progress: " << finished << "/" << total
                     << " (" << percent << "%)" << flush;
                if (finished >= total) break;
                this_thread::sleep_for(std::chrono::milliseconds(200));
            }
            cerr << '\n';
        });
    }
    for (unsigned int t = 0; t < threads; ++t) {
        pool.emplace_back(worker);
    }
    for (auto& t : pool) t.join();
    if (monitor.joinable()) monitor.join();

    u128 total = 0;
    for (u64 v : results) total += v;
    return static_cast<u64>(total);
}

bool validate(const vector<int>& primes, const vector<int>& max_exp,
              const vector<vector<vector<Point>>>& reps) {
    array<uint8_t, 8> e0{}; // all zero
    u64 f1 = compute_f(1, e0, reps);
    if (f1 != 1) {
        cerr << "Validation failed: f(1) != 1\n";
        return false;
    }
    array<uint8_t, 8> e5 = e0;
    e5[0] = 1; // 5^1
    if (compute_f(5, e5, reps) != 38) {
        cerr << "Validation failed: f(5) != 38\n";
        return false;
    }
    array<uint8_t, 8> e25 = e0;
    e25[0] = 2; // 5^2
    if (compute_f(25, e25, reps) != 167) {
        cerr << "Validation failed: f(25) != 167\n";
        return false;
    }
    array<uint8_t, 8> e13 = e0;
    e13[1] = 1; // 13^1
    if (compute_f(13, e13, reps) != 38) {
        cerr << "Validation failed: f(13) != 38\n";
        return false;
    }

    vector<int> exp_325 = {2, 1, 0, 0, 0, 0, 0, 0};
    if (compute_S(primes, exp_325, reps, false, false) != 2370) {
        cerr << "Validation failed: S(325) != 2370\n";
        return false;
    }
    vector<int> exp_1105 = {1, 1, 1, 0, 0, 0, 0, 0};
    if (compute_S(primes, exp_1105, reps, false, false) != 5535) {
        cerr << "Validation failed: S(1105) != 5535\n";
        return false;
    }
    return true;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const vector<int> primes = {5, 13, 17, 29, 37, 41, 53, 61};
    const vector<int> max_exp = {6, 3, 2, 1, 1, 1, 1, 1};

    auto reps = precompute_reps(primes, max_exp);

    if (!validate(primes, max_exp, reps)) return 1;

    u64 answer = compute_S(primes, max_exp, reps, true, true);
    cout << answer << "\n";
    return 0;
}
