#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;

static constexpr int MOD = 1'000'000'007;

static inline int mod_add(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    if (s < 0) s += MOD;
    return s;
}

static inline int mod_mul(i64 a, i64 b) { return (int)((a % MOD) * (b % MOD) % MOD); }

static int mod_pow(int a, int e) {
    long long r = 1, x = a;
    while (e > 0) {
        if (e & 1) r = (r * x) % MOD;
        x = (x * x) % MOD;
        e >>= 1;
    }
    return (int)r;
}

struct Comb {
    std::vector<int> fact, ifact;

    explicit Comb(int n) : fact(n + 1, 1), ifact(n + 1, 1) {
        for (int i = 1; i <= n; ++i) fact[i] = mod_mul(fact[i - 1], i);
        ifact[n] = mod_pow(fact[n], MOD - 2);
        for (int i = n; i >= 1; --i) ifact[i - 1] = mod_mul(ifact[i], i);
    }

    int C(int n, int k) const {
        if (k < 0 || k > n) return 0;
        return mod_mul(fact[n], mod_mul(ifact[k], ifact[n - k]));
    }
};

static int F30(int n, const Comb &cb) {
    // Hilbert series: H(t)=sum F(30,n)t^n = P(t)/(1-t)^11
    static constexpr int Pdeg = 15;
    static const int p[Pdeg + 1] = {
        1, 19, 33, 6, 0, 0, 0, 0, 0, 0, 0, -3, -33, -26, -6, 7,
    };

    int ans = 0;
    for (int k = 0; k <= Pdeg; ++k) {
        if (k > n) break;
        int ways = cb.C(n - k + 10, 10);
        ans = mod_add(ans, mod_mul(p[k], ways));
    }
    return ans;
}

static int F_bruteforce(int m, int n) {
    std::unordered_set<u64> seen;
    seen.reserve(1u << 20);
    std::vector<u64> cur{1};
    for (int step = 0; step < n; ++step) {
        std::unordered_set<u64> nxt;
        nxt.reserve(cur.size() * (size_t)m);
        for (u64 x : cur) {
            for (int a = 1; a <= m; ++a) nxt.insert(x * (u64)a);
        }
        cur.assign(nxt.begin(), nxt.end());
    }
    return (int)cur.size();
}

int main() {
    const int N = 10'001;
    Comb cb(N + 10);

    assert(F_bruteforce(9, 2) == 36);
    assert(F_bruteforce(30, 2) == 308);
    assert(F_bruteforce(30, 3) == 1909);
    assert(F_bruteforce(30, 4) == 8679);

    assert(F30(0, cb) == 1);
    assert(F30(1, cb) == 30);
    assert(F30(2, cb) == 308);
    assert(F30(3, cb) == 1909);
    assert(F30(4, cb) == 8679);

    std::cout << F30(N, cb) << "\n";
    return 0;
}

