#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using u64 = uint64_t;
using u128 = unsigned __int128;
using u256 = boost::multiprecision::uint256_t;

static inline int popcount_u64(u64 x) {
    return __builtin_popcountll(x);
}

struct Fenwick {
    int n;
    vector<int> bit; // 1-indexed

    explicit Fenwick(int n_) : n(n_), bit(n_ + 1, 0) {}

    void add(int idx0, int delta) {
        for (int i = idx0 + 1; i <= n; i += i & -i) bit[i] += delta;
    }

    int sum_prefix(int idx0) const {
        // sum over [0, idx0)
        int res = 0;
        for (int i = idx0; i > 0; i -= i & -i) res += bit[i];
        return res;
    }

    int total() const {
        return sum_prefix(n);
    }

    int kth(int k) const {
        // 0-based: smallest idx such that prefix_sum(idx+1) > k
        int idx = 0;
        int bitMask = 1;
        while ((bitMask << 1) <= n) bitMask <<= 1;
        for (int step = bitMask; step; step >>= 1) {
            int next = idx + step;
            if (next <= n && bit[next] <= k) {
                k -= bit[next];
                idx = next;
            }
        }
        return idx; // 0-based
    }
};

static inline u256 quad_score(u64 a, u64 b, u64 c, u64 d) {
    // a<=b<=c<=d, return 16*A^2 = Π(sum-2*side); 0 if invalid.
    if (d >= a + b + c) return 0;
    u64 s = a + b + c + d;
    u64 f1 = s - 2 * a;
    u64 f2 = s - 2 * b;
    u64 f3 = s - 2 * c;
    u64 f4 = s - 2 * d;
    // All factors are positive given the inequality above.
    u256 P = (u256)f1;
    P *= (u256)f2;
    P *= (u256)f3;
    P *= (u256)f4;
    return P;
}

static void print_u128(u128 x) {
    if (x == 0) {
        cout << 0;
        return;
    }
    string s;
    while (x > 0) {
        int digit = (int)(x % 10);
        s.push_back(char('0' + digit));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    cout << s;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int N = 3'000'000;

    // Precompute powers needed for u_n.
    // For n<=3e6: popcount(3n) <= 24, popcount(2n) <= 23.
    u64 pow2[25];
    for (int i = 0; i <= 24; i++) pow2[i] = 1ULL << i;
    u64 pow3[24];
    pow3[0] = 1;
    for (int i = 1; i <= 23; i++) pow3[i] = pow3[i - 1] * 3ULL;

    // Generate sequence U_1..U_N.
    vector<u64> U;
    U.reserve(N);
    for (int n = 1; n <= N; n++) {
        int b3 = popcount_u64(3ULL * (u64)n);
        int b2 = popcount_u64(2ULL * (u64)n);
        int b1 = popcount_u64((u64)n + 1ULL);
        u64 u = pow2[b3] + pow3[b2] + (u64)b1;
        U.push_back(u);
    }

    // Coordinate-compress values.
    vector<u64> vals = U;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());

    Fenwick fw((int)vals.size());
    vector<int> cnt(vals.size(), 0);

    u256 bestScore = 0;
    u64 bestPer = 0;
    u128 ans = 0;

    u128 check_sum_4_150 = 0;

    for (int n = 1; n <= N; n++) {
        u64 v = U[n - 1];
        int idx = (int)(lower_bound(vals.begin(), vals.end(), v) - vals.begin());

        int prefixLess = fw.sum_prefix(idx);
        int pos = prefixLess + cnt[idx]; // place at the end of the equal-value block

        cnt[idx]++;
        fw.add(idx, 1);

        if (n >= 4) {
            int lo = max(0, pos - 3);
            int hi = min(n - 1, pos + 3);
            vector<u64> win;
            win.reserve(hi - lo + 1);
            for (int p = lo; p <= hi; p++) {
                int id = fw.kth(p);
                win.push_back(vals[id]);
            }

            int start_min = max(0, pos - 3);
            int start_max = min(pos, n - 4);
            for (int s = start_min; s <= start_max; s++) {
                int off = s - lo;
                u64 a = win[off + 0];
                u64 b = win[off + 1];
                u64 c = win[off + 2];
                u64 d = win[off + 3];
                u256 sc = quad_score(a, b, c, d);
                if (sc == 0) continue;
                u64 per = a + b + c + d;
                if (sc > bestScore || (sc == bestScore && per > bestPer)) {
                    bestScore = sc;
                    bestPer = per;
                }
            }

            ans += (u128)bestPer;

            if (n == 5) assert(bestPer == 59);
            if (n == 10) assert(bestPer == 118);
            if (n == 150) assert(bestPer == 3223);
            if (n <= 150) check_sum_4_150 += (u128)bestPer;
        }
    }

    // sum_{4<=n<=150} f(U_n)
    assert((u64)check_sum_4_150 == 234761ULL);

    print_u128(ans);
    cout << "\n";
    return 0;
}
