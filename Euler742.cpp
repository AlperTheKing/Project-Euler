#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = std::int64_t;

struct Edge {
    int x;
    int y;
};

static bool slope_less(const Edge& a, const Edge& b) {
    return static_cast<i64>(a.x) * b.y < static_cast<i64>(b.x) * a.y;
}

static std::vector<Edge> get_edge_candidates(int num_edges) {
    std::vector<std::vector<int>> primitive(num_edges + 1, std::vector<int>(num_edges + 1, 0));
    for (int y = 1; y <= num_edges; ++y) {
        for (int x = 1; x < y; ++x) {
            if (std::gcd(x, y) == 1) primitive[y][x] = 1;
        }
    }

    std::vector<std::vector<int>> pref(num_edges + 1, std::vector<int>(num_edges + 1, 0));
    for (int y = 1; y <= num_edges; ++y) {
        int row = 0;
        for (int x = 1; x <= num_edges; ++x) {
            row += primitive[y][x];
            pref[y][x] = pref[y - 1][x] + row;
        }
    }

    std::vector<Edge> edges;
    edges.reserve(num_edges * num_edges / 2);
    for (int y = 1; y <= num_edges; ++y) {
        for (int x = 1; x < y; ++x) {
            if (!primitive[y][x]) continue;
            const int count_smaller_edges = pref[y][x];
            if (2 + count_smaller_edges > num_edges) break;
            edges.push_back({x, y});
        }
    }
    std::sort(edges.begin(), edges.end(), slope_less);
    return edges;
}

using SAPair = std::pair<int, i64>;

static std::vector<SAPair> filter_convex_hull(std::vector<SAPair> domain) {
    std::sort(domain.begin(), domain.end());
    std::vector<SAPair> out;
    out.reserve(domain.size());

    for (const auto& cur : domain) {
        const int s = cur.first;
        const i64 a = cur.second;

        if (!out.empty() && out.back().first == s) continue;

        while (out.size() >= 2) {
            const int s1 = out[out.size() - 1].first;
            const i64 a1 = out[out.size() - 1].second;
            const int s2 = out[out.size() - 2].first;
            const i64 a2 = out[out.size() - 2].second;

            __int128 lhs = static_cast<__int128>(a - a1) * static_cast<__int128>(s - s2);
            __int128 rhs = static_cast<__int128>(a - a2) * static_cast<__int128>(s - s1);
            if (lhs >= rhs) break;
            out.pop_back();
        }

        if (!out.empty() && out.back().second <= a) continue;
        out.push_back(cur);
    }
    return out;
}

static std::vector<SAPair> combine_convex_hulls(const std::vector<SAPair>& d0,
                                                const std::vector<SAPair>& d1) {
    std::vector<SAPair> merged;
    merged.reserve(d0.size() + d1.size());
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < d0.size() && j < d1.size()) {
        if (d0[i] <= d1[j]) {
            merged.push_back(d0[i++]);
        } else {
            merged.push_back(d1[j++]);
        }
    }
    while (i < d0.size()) merged.push_back(d0[i++]);
    while (j < d1.size()) merged.push_back(d1[j++]);
    return filter_convex_hull(std::move(merged));
}

static std::vector<SAPair> update_sa_domain(const std::vector<SAPair>& domain, const Edge& edge) {
    std::vector<SAPair> out;
    out.reserve(domain.size());
    const int x = edge.x;
    const int y = edge.y;
    for (const auto& [s, a] : domain) {
        const int ns = s + 2 * y;
        const i64 na = a + 2LL * x * (s + y);
        out.push_back({ns, na});
    }
    return out;
}

static i64 solve(int num_edges) {
    if (num_edges < 4 || (num_edges & 3) != 0) return -1;
    const int q_edges = num_edges / 4;
    if (q_edges == 1) return 1;

    std::vector<std::vector<SAPair>> domain(q_edges + 1);
    domain[1].push_back({1, 0});

    const auto edges = get_edge_candidates(q_edges);
    for (const auto& edge : edges) {
        for (int used = q_edges - 1; used >= 1; --used) {
            if (domain[used].empty()) continue;
            auto with_new = update_sa_domain(domain[used], edge);
            domain[used + 1] = combine_convex_hulls(domain[used + 1], with_new);
        }
    }

    i64 best = (1LL << 62);
    for (int left = 1; left < q_edges; ++left) {
        const auto& d0 = domain[left];
        const auto& d1 = domain[q_edges - left];
        for (const auto& [s0, a0] : d0) {
            for (const auto& [s1, a1] : d1) {
                const i64 val = a0 + a1 + static_cast<i64>(s0 + 2) * static_cast<i64>(s1 + 2) - 2;
                if (val < best) best = val;
            }
        }
    }
    return best;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (solve(4) != 1) {
        std::cerr << "Validation failed: A(4) != 1\n";
        return 1;
    }
    if (solve(8) != 7) {
        std::cerr << "Validation failed: A(8) != 7\n";
        return 1;
    }
    if (solve(40) != 1039) {
        std::cerr << "Validation failed: A(40) != 1039\n";
        return 1;
    }
    if (solve(100) != 17473) {
        std::cerr << "Validation failed: A(100) != 17473\n";
        return 1;
    }

    std::cout << solve(1000) << '\n';
    return 0;
}
