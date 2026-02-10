#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u128 = __uint128_t;

static constexpr int MOD = 1'000'000'007;

static inline int add_mod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

int main() {
    constexpr int N = 2000;
    const int maxDepth = (N - 1) / 5;

    std::vector<int> next(N + 1, 0), cur(N + 1, 0);

    for (int d = maxDepth; d >= 0; --d) {
        const int S = N - 5 * d;
        std::fill(cur.begin(), cur.end(), 0);

        const int maxK = S - 2;
        std::vector<u128> conv((maxK >= 0) ? (size_t)maxK + 1 : 0, 0);

        for (int s = 1; s <= S; ++s) {
            long long val = 0;
            if (s == 1) val += d;
            if (s >= 6) val += next[s - 5];
            if (s >= 4) val += (int)(conv[s - 2] % (u128)MOD);
            cur[s] = (int)(val % MOD);

            if (maxK < s + 1) continue;
            const int u_max = std::min(s, maxK - s);
            for (int u = 1; u <= u_max; ++u) {
                u128 add = (u128)cur[s] * (u128)cur[u];
                if (u != s) add *= 2;
                conv[s + u] += add;
            }
        }

        next.swap(cur);
    }

    std::vector<int> pref(N + 1, 0);
    for (int s = 1; s <= N; ++s) {
        pref[s] = add_mod(pref[s - 1], next[s]);
    }

    assert(pref[6] == 1);
    assert(pref[9] == 2);
    assert(pref[15] == 20);
    assert(pref[35] == 3166438);

    std::cout << pref[N] << "\n";
    return 0;
}

