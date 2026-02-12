#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;

static u64 highest_pow2_leq(u64 x) {
    return 1ULL << (63 - __builtin_clzll(x));
}

static u64 f_fast(u64 n, u64 k) {
    u64 cur = k;
    u64 sum = cur;
    while (cur < n) {
        u64 step = highest_pow2_leq(n - cur);
        cur += step;
        sum += cur;
    }
    return sum;
}

static u64 f_bruteforce(int n, int k) {
    std::vector<int> parent(n + 1, 0);
    std::vector<std::vector<int>> children(n + 1);
    int root = 1;

    for (int x = 2; x <= n; ++x) {
        std::vector<int> path;
        int u = root;
        while (true) {
            path.push_back(u);
            if (children[u].empty()) break;
            int best = children[u][0];
            for (int c : children[u]) best = std::max(best, c);
            u = best;
        }

        for (int i = 1; i < static_cast<int>(path.size()); ++i) {
            int p = path[i - 1], c = path[i];
            auto& v = children[p];
            v.erase(std::find(v.begin(), v.end(), c));
            parent[c] = 0;
        }

        root = x;
        for (int v : path) {
            parent[v] = x;
            children[x].push_back(v);
        }
    }

    u64 sum = 0;
    int u = k;
    while (u != 0) {
        sum += static_cast<u64>(u);
        u = parent[u];
    }
    return sum;
}

static u64 pow_u64(u64 a, int e) {
    u64 r = 1;
    while (e-- > 0) r *= a;
    return r;
}

int main() {
    assert(f_fast(6, 1) == 12ULL);
    assert(f_fast(10, 3) == 29ULL);

    for (int n = 2; n <= 25; ++n) {
        for (int k = 1; k <= n; ++k) {
            assert(f_fast(n, k) == f_bruteforce(n, k));
        }
    }

    u64 n = pow_u64(10, 17);
    u64 k = pow_u64(9, 17);
    std::cout << f_fast(n, k) << '\n';
    return 0;
}
