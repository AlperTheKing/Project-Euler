#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Dinic {
    struct Edge {
        int to;
        int rev;
        long double cap;
    };

    std::vector<std::vector<Edge>> g;
    std::vector<int> level;
    std::vector<int> it;

    explicit Dinic(int n) : g(n), level(n), it(n) {}

    void add_edge(int u, int v, long double cap) {
        Edge a{v, static_cast<int>(g[v].size()), cap};
        Edge b{u, static_cast<int>(g[u].size()), 0.0L};
        g[u].push_back(a);
        g[v].push_back(b);
    }

    bool bfs(int s, int t) {
        std::fill(level.begin(), level.end(), -1);
        std::queue<int> q;
        level[s] = 0;
        q.push(s);
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (const Edge& e : g[v]) {
                if (e.cap > 1e-18L && level[e.to] < 0) {
                    level[e.to] = level[v] + 1;
                    q.push(e.to);
                }
            }
        }
        return level[t] >= 0;
    }

    long double dfs(int v, int t, long double f) {
        if (v == t) return f;
        for (int& i = it[v]; i < static_cast<int>(g[v].size()); ++i) {
            Edge& e = g[v][i];
            if (e.cap <= 1e-18L || level[e.to] != level[v] + 1) continue;
            long double got = dfs(e.to, t, std::min(f, e.cap));
            if (got > 1e-18L) {
                e.cap -= got;
                g[e.to][e.rev].cap += got;
                return got;
            }
        }
        return 0.0L;
    }

    long double max_flow(int s, int t) {
        long double flow = 0.0L;
        while (bfs(s, t)) {
            std::fill(it.begin(), it.end(), 0);
            while (true) {
                long double f = dfs(s, t, std::numeric_limits<long double>::infinity());
                if (f <= 1e-18L) break;
                flow += f;
            }
        }
        return flow;
    }
};

static std::vector<int> sieve_primes(int n) {
    std::vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) spf[i] = i;
    for (int i = 2; static_cast<u64>(i) * i <= static_cast<u64>(n); ++i) {
        if (spf[i] != i) continue;
        for (int j = i * i; j <= n; j += i) {
            if (spf[j] == j) spf[j] = i;
        }
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (spf[i] == i) primes.push_back(i);
    }
    return primes;
}

static long double solve_ln(int n) {
    const std::vector<int> primes = sieve_primes(n);

    std::vector<int> a7;
    std::vector<int> b9;
    long double mandatory_log = 0.0L;
    for (int p : primes) {
        const int r = p % 10;
        if (r == 3) mandatory_log += std::log(static_cast<long double>(p));
        if (r == 7) a7.push_back(p);
        if (r == 9) b9.push_back(p);
    }

    std::vector<bool> forced(n + 1, false);
    long double forced_log = 0.0L;
    for (int p : a7) {
        u128 cube = static_cast<u128>(p) * p * p;
        if (cube <= static_cast<u128>(n)) {
            forced[p] = true;
            forced_log += std::log(static_cast<long double>(p));
        }
    }

    std::vector<int> nonforced_a;
    nonforced_a.reserve(a7.size());
    for (int p : a7) {
        if (!forced[p]) nonforced_a.push_back(p);
    }

    std::vector<std::pair<int, int>> edges;
    edges.reserve(10000);
    std::vector<bool> used_a(n + 1, false);
    std::vector<bool> used_b(n + 1, false);

    for (int a : nonforced_a) {
        int limit = n / a;
        auto it = std::upper_bound(b9.begin(), b9.end(), limit);
        for (auto jt = b9.begin(); jt != it; ++jt) {
            int b = *jt;
            edges.push_back({a, b});
            used_a[a] = true;
            used_b[b] = true;
        }
    }

    long double cover_log = 0.0L;
    if (!edges.empty()) {
        std::vector<int> id_a(n + 1, -1);
        std::vector<int> id_b(n + 1, -1);

        int node_count = 1;
        for (int a : nonforced_a) {
            if (used_a[a]) id_a[a] = node_count++;
        }
        for (int b : b9) {
            if (used_b[b]) id_b[b] = node_count++;
        }
        const int source = 0;
        const int sink = node_count++;

        Dinic dinic(node_count);
        const long double inf_cap = 1e30L;

        for (int a : nonforced_a) {
            if (used_a[a]) {
                dinic.add_edge(source, id_a[a], std::log(static_cast<long double>(a)));
            }
        }
        for (int b : b9) {
            if (used_b[b]) {
                dinic.add_edge(id_b[b], sink, std::log(static_cast<long double>(b)));
            }
        }
        for (const auto& e : edges) {
            dinic.add_edge(id_a[e.first], id_b[e.second], inf_cap);
        }

        cover_log = dinic.max_flow(source, sink);
    }

    return mandatory_log + forced_log + cover_log;
}

static long double round6(long double x) {
    return std::round(x * 1'000'000.0L) / 1'000'000.0L;
}

int main() {
    assert(std::fabsl(round6(solve_ln(40)) - 6.799056L) < 1e-12L);
    assert(std::fabsl(round6(solve_ln(2800)) - 715.019337L) < 1e-12L);

    std::cout << std::fixed << std::setprecision(6) << round6(solve_ln(1'000'000)) << '\n';
    return 0;
}
