#include <algorithm>
#include <cassert>
#include <cstdint>
#include <deque>
#include <iostream>
#include <vector>

using i64 = std::int64_t;

static int path_max(const std::vector<int>& w, int l, int r) {
    if (l > r) return 0;
    int a = 0, b = 0;
    for (int i = l; i <= r; ++i) {
        int c = std::max(b, a + w[i]);
        a = b;
        b = c;
    }
    return b;
}

static int D_value(int n) {
    std::vector<int> to(n), indeg(n, 0);
    std::vector<std::vector<int>> preds(n);

    for (int x = 0; x < n; ++x) {
        i64 y = x;
        int v = static_cast<int>(((y * y % n) * y + y + 1) % n);
        to[x] = v;
        ++indeg[v];
        preds[v].push_back(x);
    }

    std::deque<int> q;
    for (int i = 0; i < n; ++i) {
        if (indeg[i] == 0) q.push_back(i);
    }

    std::vector<int> order;
    order.reserve(n);
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        order.push_back(u);

        int v = to[u];
        --indeg[v];
        if (indeg[v] == 0) q.push_back(v);
    }

    std::vector<int> dp0(n, 0), dp1(n, 0);
    for (int t = 0; t < static_cast<int>(order.size()); ++t) {
        int u = order[t];
        int base = 0;
        int best = 0;

        for (int c : preds[u]) {
            if (indeg[c] == 0) {
                base += dp0[c];
                best = std::max(best, dp1[c] - dp0[c]);
            }
        }

        dp0[u] = base + std::max(0, best);
        dp1[u] = base + 1;
    }

    std::vector<char> vis(n, 0);
    int ans = 0;

    for (int s = 0; s < n; ++s) {
        if (indeg[s] == 0 || vis[s]) continue;

        std::vector<int> cyc;
        int u = s;
        while (!vis[u]) {
            vis[u] = 1;
            cyc.push_back(u);
            u = to[u];
        }

        int k = static_cast<int>(cyc.size());
        std::vector<int> base(k, 0), gain(k, 0);

        for (int i = 0; i < k; ++i) {
            int x = cyc[i];
            int b = 0;
            int g = 0;

            for (int c : preds[x]) {
                if (indeg[c] == 0) {
                    b += dp0[c];
                    g = std::max(g, dp1[c] - dp0[c]);
                }
            }

            base[i] = b;
            gain[i] = std::max(0, g);
        }

        int comp = 0;
        for (int i = 0; i < k; ++i) comp += base[i] + gain[i];

        if (k == 1) {
            ans += comp;
            continue;
        }

        std::vector<int> w(k, 0);
        for (int i = 0; i < k; ++i) {
            int j = (i + 1 == k ? 0 : i + 1);
            w[i] = 1 - gain[i] - gain[j];
        }

        int bestA = path_max(w, 1, k - 1);
        int bestB = w[0] + path_max(w, 2, k - 2);
        int best = std::max(0, std::max(bestA, bestB));

        ans += comp + best;
    }

    return ans;
}

int main() {
    assert(D_value(5) == 1);
    assert(D_value(10) == 3);

    int sum = 0;
    for (int i = 1; i <= 100; ++i) {
        sum += D_value(100'000 + i);
    }

    std::cout << sum << '\n';
    return 0;
}
