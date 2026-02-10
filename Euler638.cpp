#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;

static constexpr int MOD = 1'000'000'007;

static inline int mod_add(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    if (s < 0) s += MOD;
    return s;
}

static inline int mod_mul(i64 a, i64 b) { return (int)(a * b % MOD); }

static int mod_pow(int a, i64 e) {
    i64 r = 1, x = a;
    while (e > 0) {
        if (e & 1) r = (r * x) % MOD;
        x = (x * x) % MOD;
        e >>= 1;
    }
    return (int)r;
}

static int mod_inv(int a) { return mod_pow(a, MOD - 2); }

static int binom_small(int n, int k) {
    static std::vector<int> fact, ifact;
    if (fact.empty()) {
        const int N = 50;
        fact.assign(N + 1, 1);
        for (int i = 1; i <= N; ++i) fact[i] = mod_mul(fact[i - 1], i);
        ifact.assign(N + 1, 1);
        ifact[N] = mod_inv(fact[N]);
        for (int i = N; i >= 1; --i) ifact[i - 1] = mod_mul(ifact[i], i);
    }
    if (k < 0 || k > n) return 0;
    return mod_mul(fact[n], mod_mul(ifact[k], ifact[n - k]));
}

static int C(int a, int b, int q) {
    if (a < 0 || b < 0) return 0;
    if (a == 0 || b == 0) return 1;
    if (q == 1) return binom_small(a + b, a);

    int k = a < b ? a : b;
    int m = a + b - k;
    int qm = mod_pow(q, m);

    int denom = 1;
    int numer = 1;
    int qj = 1;
    for (int j = 1; j <= k; ++j) {
        qj = mod_mul(qj, q);  // q^j
        denom = mod_mul(denom, mod_add(1, -qj));
        numer = mod_mul(numer, mod_add(1, -mod_mul(qm, qj)));
    }

    assert(denom != 0);
    return mod_mul(numer, mod_inv(denom));
}

int main() {
    assert(C(2, 2, 1) == 6);
    assert(C(2, 2, 2) == 35);
    assert(C(10, 10, 1) == 184'756);
    assert(C(15, 10, 3) == 880'419'838);
    assert(C(10'000, 10'000, 4) == 395'913'804);

    i64 ans = 0;
    i64 ten = 1;
    for (int k = 1; k <= 7; ++k) {
        ten *= 10;
        const int n = (int)(ten + k);
        ans = mod_add((int)ans, C(n, n, k));
    }
    std::cout << (int)ans << "\n";
    return 0;
}

