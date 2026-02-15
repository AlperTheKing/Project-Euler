#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct State {
    std::array<std::uint8_t, 7> labels{};
    std::array<std::uint8_t, 7> sizes{};
    std::uint8_t k = 0;
    std::uint8_t best = 0;

    bool operator==(const State& other) const {
        return k == other.k && best == other.best && labels == other.labels && sizes == other.sizes;
    }
};

struct StateHash {
    std::size_t operator()(const State& s) const noexcept {
        std::uint64_t h = 1469598103934665603ULL;
        const auto mix = [&](const std::uint8_t x) {
            h ^= static_cast<std::uint64_t>(x);
            h *= 1099511628211ULL;
        };

        for (const std::uint8_t x : s.labels) {
            mix(x);
        }
        for (const std::uint8_t x : s.sizes) {
            mix(x);
        }
        mix(s.k);
        mix(s.best);
        return static_cast<std::size_t>(h);
    }
};

struct DSU {
    std::array<int, 14> parent{};
    std::array<int, 14> comp_size{};

    explicit DSU(const int n) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
            comp_size[i] = 0;
        }
    }

    int find(const int x) {
        int r = x;
        while (parent[r] != r) {
            r = parent[r];
        }
        int y = x;
        while (parent[y] != y) {
            const int p = parent[y];
            parent[y] = r;
            y = p;
        }
        return r;
    }

    void unite(const int a, const int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) {
            return;
        }
        if (comp_size[ra] < comp_size[rb]) {
            std::swap(ra, rb);
        }
        parent[rb] = ra;
        comp_size[ra] += comp_size[rb];
    }
};

State transition(const State& s, const int w, const int mask) {
    const int old_k = static_cast<int>(s.k);
    std::array<int, 7> cell_node{};
    cell_node.fill(-1);

    int new_cells = 0;
    for (int c = 0; c < w; ++c) {
        if ((mask >> c) & 1) {
            cell_node[static_cast<std::size_t>(c)] = old_k + new_cells;
            ++new_cells;
        }
    }

    DSU dsu(old_k + new_cells);

    for (int i = 0; i < old_k; ++i) {
        dsu.comp_size[static_cast<std::size_t>(i)] = static_cast<int>(s.sizes[static_cast<std::size_t>(i)]);
    }
    for (int c = 0; c < w; ++c) {
        if (cell_node[static_cast<std::size_t>(c)] != -1) {
            dsu.comp_size[static_cast<std::size_t>(cell_node[static_cast<std::size_t>(c)])] = 1;
        }
    }

    for (int c = 0; c < w; ++c) {
        const int node = cell_node[static_cast<std::size_t>(c)];
        if (node == -1) {
            continue;
        }

        const int above_label = static_cast<int>(s.labels[static_cast<std::size_t>(c)]);
        if (above_label > 0) {
            dsu.unite(node, above_label - 1);
        }
        if (c > 0) {
            const int left_node = cell_node[static_cast<std::size_t>(c - 1)];
            if (left_node != -1) {
                dsu.unite(node, left_node);
            }
        }
    }

    State out;
    out.best = s.best;

    std::array<int, 14> root_to_new{};
    root_to_new.fill(0);
    int next_label = 0;

    for (int c = 0; c < w; ++c) {
        const int node = cell_node[static_cast<std::size_t>(c)];
        if (node == -1) {
            out.labels[static_cast<std::size_t>(c)] = 0;
            continue;
        }

        const int root = dsu.find(node);
        int nl = root_to_new[static_cast<std::size_t>(root)];
        if (nl == 0) {
            ++next_label;
            nl = next_label;
            root_to_new[static_cast<std::size_t>(root)] = nl;
            const int sz = dsu.comp_size[static_cast<std::size_t>(root)];
            out.sizes[static_cast<std::size_t>(nl - 1)] = static_cast<std::uint8_t>(sz);
            if (sz > static_cast<int>(out.best)) {
                out.best = static_cast<std::uint8_t>(sz);
            }
        }
        out.labels[static_cast<std::size_t>(c)] = static_cast<std::uint8_t>(nl);
    }

    out.k = static_cast<std::uint8_t>(next_label);
    return out;
}

long double expected_max_area(const int w, const int h) {
    std::unordered_map<State, u64, StateHash> cur;
    std::unordered_map<State, u64, StateHash> next;
    cur.reserve(4096);
    next.reserve(4096);

    cur[State{}] = 1ULL;

    const int row_masks = 1 << w;
    for (int r = 0; r < h; ++r) {
        next.clear();
        next.reserve(static_cast<std::size_t>(cur.size()) << w);
        for (const auto& kv : cur) {
            const State& s = kv.first;
            const u64 ways = kv.second;

            for (int mask = 0; mask < row_masks; ++mask) {
                const State ns = transition(s, w, mask);
                next[ns] += ways;
            }
        }
        cur.swap(next);
        next.reserve(std::max<std::size_t>(cur.size() << w, 4096ULL));
    }

    u128 numerator = 0;
    for (const auto& kv : cur) {
        numerator += static_cast<u128>(kv.first.best) * static_cast<u128>(kv.second);
    }

    const long double denominator = std::ldexp(1.0L, w * h);
    return static_cast<long double>(numerator) / denominator;
}

}  // namespace

int main() {
    const long double e22 = expected_max_area(2, 2);
    assert(std::fabsl(e22 - 1.875L) < 1e-12L);

    const long double e44 = expected_max_area(4, 4);
    assert(std::fabsl(e44 - 5.76487732L) < 5e-9L);

    const long double e77 = expected_max_area(7, 7);
    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(e77) << '\n';
    return 0;
}
