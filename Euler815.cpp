#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>

using i64 = std::int64_t;

static constexpr int kBase = 61;

static inline int pack(int n0, int n1, int n2, int n3) {
    return (((n0 * kBase) + n1) * kBase + n2) * kBase + n3;
}

static inline void unpack(int key, int& n0, int& n1, int& n2, int& n3) {
    n3 = key % kBase;
    key /= kBase;
    n2 = key % kBase;
    key /= kBase;
    n1 = key % kBase;
    key /= kBase;
    n0 = key;
}

static long double prob_max_lt(int n, int k) {
    std::unordered_map<int, long double> cur;
    cur.reserve(200000);
    cur[pack(n, 0, 0, 0)] = 1.0L;

    const int total_cards = 4 * n;

    for (int step = 0; step < total_cards; ++step) {
        std::unordered_map<int, long double> nxt;
        nxt.reserve(cur.size() * 2 + 64);

        for (const auto& it : cur) {
            const int key = it.first;
            const long double prob = it.second;

            int n0, n1, n2, n3;
            unpack(key, n0, n1, n2, n3);

            const int remaining = 4 * n0 + 3 * n1 + 2 * n2 + n3;
            const int piles = n1 + n2 + n3;
            const long double inv = 1.0L / static_cast<long double>(remaining);

            if (n0 > 0 && piles + 1 < k) {
                const long double add = prob * static_cast<long double>(4 * n0) * inv;
                nxt[pack(n0 - 1, n1 + 1, n2, n3)] += add;
            }
            if (n1 > 0 && piles < k) {
                const long double add = prob * static_cast<long double>(3 * n1) * inv;
                nxt[pack(n0, n1 - 1, n2 + 1, n3)] += add;
            }
            if (n2 > 0 && piles < k) {
                const long double add = prob * static_cast<long double>(2 * n2) * inv;
                nxt[pack(n0, n1, n2 - 1, n3 + 1)] += add;
            }
            if (n3 > 0 && piles - 1 < k) {
                const long double add = prob * static_cast<long double>(n3) * inv;
                nxt[pack(n0, n1, n2, n3 - 1)] += add;
            }
        }

        cur.swap(nxt);
        if (cur.empty()) {
            return 0.0L;
        }
    }

    auto it = cur.find(pack(0, 0, 0, 0));
    return (it == cur.end() ? 0.0L : it->second);
}

static long double expected_max_piles(int n) {
    long double ans = 0.0L;
    for (int k = 1; k <= n; ++k) {
        ans += (1.0L - prob_max_lt(n, k));
    }
    return ans;
}

int main() {
    const long double e2 = expected_max_piles(2);
    assert(static_cast<i64>(std::llround(e2 * 100000000.0L)) == 197142857LL);

    const long double e60 = expected_max_piles(60);
    const long double rounded = std::round(e60 * 100000000.0L) / 100000000.0L;
    std::cout << std::fixed << std::setprecision(8) << rounded << '\n';
    return 0;
}
