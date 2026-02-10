#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

static inline u64 isqrt_u128(u128 x) {
    long double d = std::sqrt((long double)x);
    u64 r = (u64)d;
    while ((u128)(r + 1) * (r + 1) <= x) ++r;
    while ((u128)r * r > x) --r;
    return r;
}

static inline u64 isqrt_u64(u64 x) {
    u64 r = (u64)std::sqrt((long double)x);
    while ((u128)(r + 1) * (r + 1) <= x) ++r;
    while ((u128)r * r > x) --r;
    return r;
}

static inline u128 ipow_u128(u64 a, int e) {
    u128 r = 1, x = a;
    while (e > 0) {
        if (e & 1) r *= x;
        x *= x;
        e >>= 1;
    }
    return r;
}

static u64 iroot6_u128(u128 n) {
    u64 r = (u64)std::pow((long double)n, 1.0L / 6.0L);
    while (ipow_u128(r + 1, 6) <= n) ++r;
    while (ipow_u128(r, 6) > n) --r;
    return r;
}

struct Mobius {
    int N;
    std::vector<int> mu;
    std::vector<i64> pref_mu;
    std::vector<i64> pref_sqfree;

    explicit Mobius(int n) : N(n), mu(n + 1, 0), pref_mu(n + 1, 0), pref_sqfree(n + 1, 0) {
        std::vector<int> primes;
        std::vector<int> lp(n + 1, 0);
        mu[1] = 1;
        for (int i = 2; i <= n; ++i) {
            if (lp[i] == 0) {
                lp[i] = i;
                primes.push_back(i);
                mu[i] = -1;
            }
            for (int p : primes) {
                if (p > lp[i] || (i64)i * p > n) break;
                lp[i * p] = p;
                if (p == lp[i]) {
                    mu[i * p] = 0;
                    break;
                }
                mu[i * p] = -mu[i];
            }
        }
        for (int i = 1; i <= n; ++i) {
            pref_mu[i] = pref_mu[i - 1] + mu[i];
            pref_sqfree[i] = pref_sqfree[i - 1] + (mu[i] != 0);
        }
    }
};

struct MertensMemo {
    const Mobius &mb;
    std::unordered_map<u64, i64> memo;

    explicit MertensMemo(const Mobius &m) : mb(m) { memo.reserve(1 << 20); }

    i64 M(u64 n) {
        if (n <= (u64)mb.N) return mb.pref_mu[(size_t)n];
        auto it = memo.find(n);
        if (it != memo.end()) return it->second;
        i64 res = 1;
        for (u64 l = 2; l <= n;) {
            const u64 q = n / l;
            const u64 r = n / q;
            res -= (i64)(r - l + 1) * M(q);
            l = r + 1;
        }
        memo.emplace(n, res);
        return res;
    }
};

static i64 squarefree_count(u64 n, const Mobius &mb, MertensMemo &) {
    if (n <= (u64)mb.N) return mb.pref_sqfree[(size_t)n];
    const u64 r = isqrt_u64(n);
    i64 s = 0;
    for (u64 d = 1; d <= r; ++d) s += (i64)mb.mu[(size_t)d] * (i64)(n / (d * d));
    return s;
}

static i64 squarefree_even_count(u64 n, const Mobius &mb, MertensMemo &mm) {
    const i64 Q = squarefree_count(n, mb, mm);
    const i64 M = mm.M(n);
    return (Q + M) / 2;
}

static i64 f(u128 N) {
    static constexpr int LIM = 1'000'000;  // (1e9)^(2/3)
    const Mobius mb(LIM);
    MertensMemo mm(mb);

    const u64 amax = iroot6_u128(N);
    i64 ans = 0;
    for (u64 a = 1; a <= amax; ++a) {
        const u128 a6 = ipow_u128(a, 6);
        const u128 t = N / a6;
        const u64 s = isqrt_u128(t);
        const u64 bmax = isqrt_u64(s);
        ans += squarefree_even_count(bmax, mb, mm);
    }
    return ans;
}

int main() {
    {
        u128 N = 100;
        assert(f(N) == 2);
    }
    {
        u128 N = 100'000'000ULL;
        assert(f(N) == 69);
    }

    u128 N = 1;
    for (int i = 0; i < 36; ++i) N *= 10;
    std::cout << f(N) << "\n";
    return 0;
}

