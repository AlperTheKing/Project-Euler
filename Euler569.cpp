#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

// 128-bit safe multiply-add for cross products / slope comparisons.
using i128 = __int128_t;

// Upper hull maintenance for points with strictly increasing x.
static inline i128 cross(const std::vector<i64>& x, const std::vector<i64>& y, int a, int b, int c) {
    // cross( (b-a), (c-b) )
    const i128 abx = (i128)x[b] - x[a];
    const i128 aby = (i128)y[b] - y[a];
    const i128 bcx = (i128)x[c] - x[b];
    const i128 bcy = (i128)y[c] - y[b];
    return abx * bcy - aby * bcx;
}

// Compare slopes to a fixed right endpoint k:
// slope(k, j1) < slope(k, j2)  <=>  (y_k - y_{j1})/(x_k - x_{j1}) < (y_k - y_{j2})/(x_k - x_{j2})
static inline bool slope_less(const std::vector<i64>& x,
                              const std::vector<i64>& y,
                              int k,
                              int j1,
                              int j2) {
    const i128 dy1 = (i128)y[k] - y[j1];
    const i128 dx1 = (i128)x[k] - x[j1];
    const i128 dy2 = (i128)y[k] - y[j2];
    const i128 dx2 = (i128)x[k] - x[j2];
    return dy1 * dx2 < dy2 * dx1;
}

static std::vector<int> first_primes(std::size_t need) {
    // Sieve up to 1e8: pi(1e8)=5761455 > 5e6, enough for this problem.
    const int limit = 100'000'000;
    const int half = limit / 2;  // odds: n = 2*i+1 for i in [1..half]

    std::vector<bool> comp(half + 1, false);
    std::vector<int> primes;
    primes.reserve(need);
    primes.push_back(2);

    const int sqrt_limit = 10'000;  // floor(sqrt(1e8))
    for (int i = 1; (2 * i + 1) <= sqrt_limit; ++i) {
        if (comp[i]) continue;
        const int p = 2 * i + 1;
        const int start = (p * p) / 2;
        for (int j = start; j <= half; j += p) comp[j] = true;
    }

    for (int i = 1; i <= half && primes.size() < need; ++i) {
        if (!comp[i]) primes.push_back(2 * i + 1);
    }
    assert(primes.size() >= need);
    primes.resize(need);
    return primes;
}

struct Data {
    std::vector<i64> x;
    std::vector<i64> y;
    std::vector<int> hull_prev;  // previous vertex on upper hull of prefix i (valid for all i).
};

static Data build_points_and_prefix_hulls(int N) {
    // Need primes p_1..p_{2N}.
    const std::size_t need = (std::size_t)2 * (std::size_t)N;
    const std::vector<int> p = first_primes(need);

    Data d;
    d.x.assign(N + 1, 0);
    d.y.assign(N + 1, 0);
    d.hull_prev.assign(N + 1, 0);

    i64 valley_y = 0;
    i64 cur_x = 0;
    for (int k = 1; k <= N; ++k) {
        const int up = p[(std::size_t)2 * (std::size_t)k - 2];    // p_{2k-1}
        const int down = p[(std::size_t)2 * (std::size_t)k - 1];  // p_{2k}
        cur_x += up;
        d.x[k] = cur_x;
        d.y[k] = valley_y + up;
        cur_x += down;
        valley_y += (i64)up - (i64)down;
    }

    // Build upper hull for each prefix using the standard monotone chain update, but store predecessor pointers
    // so each prefix hull can be traversed by following hull_prev pointers from its rightmost vertex.
    std::vector<int> st;
    st.reserve(N);
    for (int i = 1; i <= N; ++i) {
        while (st.size() >= 2U) {
            const int a = st[st.size() - 2];
            const int b = st[st.size() - 1];
            if (cross(d.x, d.y, a, b, i) >= 0) {
                st.pop_back();
            } else {
                break;
            }
        }
        d.hull_prev[i] = st.empty() ? 0 : st.back();
        st.push_back(i);
    }

    return d;
}

static int P_full_scan(const Data& d, int k) {
    if (k <= 1) return 0;
    int cnt = 1;
    int best = k - 1;
    for (int j = k - 2; j >= 1; --j) {
        if (slope_less(d.x, d.y, k, j, best)) {
            best = j;
            ++cnt;
        }
    }
    return cnt;
}

static int P_fast(const Data& d, int k) {
    if (k <= 1) return 0;
    int cnt = 1;
    int best = k - 1;
    int t = best;

    while (t > 1) {
        // Existence test: min slope among points in [1..t-1] is achieved on the upper hull of prefix (t-1).
        bool exists = false;
        for (int v = t - 1; v != 0; v = d.hull_prev[v]) {
            if (slope_less(d.x, d.y, k, v, best)) {
                exists = true;
                break;
            }
        }
        if (!exists) break;

        int j = t - 1;
        while (j >= 1 && !slope_less(d.x, d.y, k, j, best)) --j;
        if (j <= 0) break;  // should not happen, but keeps us safe if equality ever occurs.

        ++cnt;
        best = j;
        t = j;
    }
    return cnt;
}

static u64 sum_P_range(const Data& d, int L, int R) {
    u64 acc = 0;
    for (int k = L; k <= R; ++k) acc += (u64)P_fast(d, k);
    return acc;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    constexpr int N = 2'500'000;
    const Data d = build_points_and_prefix_hulls(N);

    // Validation by direct scan for the first few peaks.
    u64 sumP = 0;
    for (int k = 1; k <= 200; ++k) {
        const int pk_fast = P_fast(d, k);
        const int pk_full = P_full_scan(d, k);
        assert(pk_fast == pk_full);
        sumP += (u64)pk_fast;
        if (k == 3) assert(pk_fast == 1);
        if (k == 9) assert(pk_fast == 3);
        if (k == 100) assert(sumP == 227);
    }

    const int start = 201;
    if (start <= N) {
        const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
        const unsigned threads = std::min<unsigned>(hw, 8u);
        std::vector<std::thread> pool;
        std::vector<u64> partial(threads, 0);

        const int total = N - start + 1;
        const int chunk = (total + (int)threads - 1) / (int)threads;
        for (unsigned ti = 0; ti < threads; ++ti) {
            const int L = start + (int)ti * chunk;
            const int R = std::min(N, L + chunk - 1);
            if (L > R) continue;
            pool.emplace_back([&d, &partial, ti, L, R]() { partial[ti] = sum_P_range(d, L, R); });
        }
        for (auto& th : pool) th.join();
        for (u64 v : partial) sumP += v;
    }

    std::cout << sumP << '\n';
    return 0;
}

