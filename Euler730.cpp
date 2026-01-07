#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

using int64 = long long;

struct Triple {
    int64 p;
    int64 q;
    int64 r;
};

static bool is_square(int64 v, int64 &root) {
    if (v < 0) return false;
    int64 r = static_cast<int64>(std::sqrt(static_cast<long double>(v)) + 0.5L);
    if (r * r == v) {
        root = r;
        return true;
    }
    return false;
}

static bool has_parent(const Triple &t) {
    const int64 p = t.p;
    const int64 q = t.q;
    const int64 r = t.r;
    const int64 x3 = -2 * p - 2 * q + 3 * r;
    if (x3 <= 0) return false;
    int64 x1 = p + 2 * q - 2 * r;
    int64 x2 = -2 * p - q + 2 * r;
    if (x1 > 0 && x2 > 0) return true;  // A^{-1}
    x2 = 2 * p + q - 2 * r;
    if (x1 > 0 && x2 > 0) return true;  // B^{-1}
    x1 = -p - 2 * q + 2 * r;
    if (x1 > 0 && x2 > 0) return true;  // C^{-1}
    return false;
}

static std::vector<std::vector<Triple>> build_roots(int m) {
    // Roots are primitive solutions with no parent under A/B/C inverses.
    // For fixed k there are finitely many roots, and p,q <= 4k suffices for k>0.
    int limit = std::max(4 * m, 10);
    std::vector<std::vector<Triple>> roots(m + 1);
    for (int k = 0; k <= m; ++k) {
        for (int p = 1; p <= limit; ++p) {
            for (int q = p; q <= limit; ++q) {
                int64 r;
                int64 v = 1LL * p * p + 1LL * q * q + k;
                if (!is_square(v, r)) continue;
                if (std::gcd<int64>(p, std::gcd<int64>(q, r)) != 1) continue;
                Triple t{p, q, r};
                if (!has_parent(t)) roots[k].push_back(t);
            }
        }
    }
    return roots;
}

struct CountPair {
    int64 total;
    int64 diag;
};

static CountPair count_from_root(const Triple &root, int64 n, std::vector<Triple> &stack) {
    stack.clear();
    stack.push_back(root);
    int64 count = 0;
    int64 diag = 0;
    while (!stack.empty()) {
        Triple cur = stack.back();
        stack.pop_back();
        int64 per = cur.p + cur.q + cur.r;
        if (per > n) continue;
        ++count;
        if (cur.p == cur.q) ++diag;

        Triple a{cur.p - 2 * cur.q + 2 * cur.r,
                 2 * cur.p - cur.q + 2 * cur.r,
                 2 * cur.p - 2 * cur.q + 3 * cur.r};
        if (a.p + a.q + a.r <= n) stack.push_back(a);

        Triple b{cur.p + 2 * cur.q + 2 * cur.r,
                 2 * cur.p + cur.q + 2 * cur.r,
                 2 * cur.p + 2 * cur.q + 3 * cur.r};
        if (b.p + b.q + b.r <= n) stack.push_back(b);

        Triple c{-cur.p + 2 * cur.q + 2 * cur.r,
                 -2 * cur.p + cur.q + 2 * cur.r,
                 -2 * cur.p + 2 * cur.q + 3 * cur.r};
        if (c.p + c.q + c.r <= n) stack.push_back(c);
    }
    return {count, diag};
}

static std::vector<int64> count_all(int m, int64 n,
                                    const std::vector<std::vector<Triple>> &roots,
                                    int threads) {
    struct Task {
        int k;
        Triple root;
    };

    std::vector<Task> tasks;
    for (int k = 0; k <= m; ++k) {
        for (const auto &r : roots[k]) tasks.push_back({k, r});
    }

    if (tasks.empty()) return std::vector<int64>(m + 1, 0);

    unsigned hw = std::thread::hardware_concurrency();
    int T = (threads > 0) ? threads : (hw ? static_cast<int>(hw) : 1);
    T = std::max(1, std::min(T, 16));
    T = std::min<int>(T, tasks.size());

    std::vector<std::vector<int64>> local(T, std::vector<int64>(m + 1, 0));
    std::atomic<size_t> next(0);

    auto worker = [&](int tid) {
        std::vector<Triple> stack;
        stack.reserve(1024);
        while (true) {
            size_t idx = next.fetch_add(1);
            if (idx >= tasks.size()) break;
            const Task &task = tasks[idx];
            int64 per = task.root.p + task.root.q + task.root.r;
            if (per > n) continue;
            CountPair cnt = count_from_root(task.root, n, stack);
            // If this component hits the diagonal p==q, swapping p and q keeps the
            // component and pairs off-diagonal nodes; adjust to unordered counts.
            if (cnt.diag > 0) {
                local[tid][task.k] += (cnt.total + cnt.diag) / 2;
            } else {
                local[tid][task.k] += cnt.total;
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(T);
    for (int t = 0; t < T; ++t) pool.emplace_back(worker, t);
    for (auto &th : pool) th.join();

    std::vector<int64> counts(m + 1, 0);
    for (int t = 0; t < T; ++t) {
        for (int k = 0; k <= m; ++k) counts[k] += local[t][k];
    }
    return counts;
}

static bool run_validation(const std::vector<std::vector<Triple>> &roots, int threads) {
    const int64 n = 10000;
    const int m = 20;
    auto counts = count_all(m, n, roots, threads);

    if (counts[0] != 703) {
        std::cerr << "Validation failed: P_0(1e4) expected 703, got " << counts[0] << "\n";
        return false;
    }
    if (counts[20] != 1979) {
        std::cerr << "Validation failed: P_20(1e4) expected 1979, got " << counts[20] << "\n";
        return false;
    }
    int64 s10 = 0;
    for (int k = 0; k <= 10; ++k) s10 += counts[k];
    if (s10 != 10956) {
        std::cerr << "Validation failed: S(10,1e4) expected 10956, got " << s10 << "\n";
        return false;
    }
    std::cerr << "Validation checkpoints passed.\n";
    return true;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const int m = 100;
    const int64 n = 100000000LL;

    unsigned hw = std::thread::hardware_concurrency();
    int threads = hw ? static_cast<int>(hw) : 1;
    threads = std::max(1, std::min(threads, 16));

    auto roots = build_roots(m);
    if (!run_validation(roots, threads)) return 1;

    auto counts = count_all(m, n, roots, threads);
    int64 total = 0;
    for (int k = 0; k <= m; ++k) total += counts[k];
    std::cout << total << "\n";
    return 0;
}
