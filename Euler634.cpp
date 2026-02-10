#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

static inline u64 isqrt_u64(u64 x) {
    u64 r = (u64)std::sqrt((long double)x);
    while ((r + 1) > 0 && (r + 1) * (u128)(r + 1) <= x) ++r;
    while (r * (u128)r > x) --r;
    return r;
}

static inline u64 icbrt_u64(u64 x) {
    u64 r = (u64)std::cbrt((long double)x);
    while ((r + 1) * (u128)(r + 1) * (u128)(r + 1) <= x) ++r;
    while (r * (u128)r * (u128)r > x) --r;
    return r;
}

static inline u64 iroot6_u64(u64 x) {
    u64 r = (u64)std::pow((long double)x, 1.0L / 6.0L);
    auto p6 = [](u64 a) -> u128 {
        u128 t = (u128)a * a;
        t *= t;          // a^4
        t *= (u128)a * a; // a^6
        return t;
    };
    while (p6(r + 1) <= x) ++r;
    while (p6(r) > x) --r;
    return r;
}

static std::vector<int> primes_upto(int n) {
    if (n < 2) return {};
    std::vector<bool> is_comp(n + 1, false);
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (!is_comp[i]) {
            primes.push_back(i);
            if ((i64)i * i <= n) {
                for (int j = i * i; j <= n; j += i) is_comp[j] = true;
            }
        }
    }
    return primes;
}

static std::vector<int> mobius_upto(int n) {
    std::vector<int> mu(n + 1, 0), primes, lp(n + 1, 0);
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
    return mu;
}

static i64 count_cubefree_upto(u64 n) {
    const int D = (int)icbrt_u64(n);
    const auto mu = mobius_upto(D);
    i64 s = 0;
    for (int d = 1; d <= D; ++d) {
        if (mu[d] == 0) continue;
        const u64 d3 = (u64)d * d * d;
        s += (i64)mu[d] * (i64)(n / d3);
    }
    return s;
}

static i64 F(u64 N) {
    const u64 vmax = icbrt_u64(N);
    std::vector<uint8_t> sqfree(vmax + 1, 1);
    sqfree[0] = 0;

    const int P = (int)isqrt_u64(vmax);
    const auto primes = primes_upto(P);
    for (int p : primes) {
        const u64 p2 = (u64)p * p;
        for (u64 m = p2; m <= vmax; m += p2) sqfree[(size_t)m] = 0;
    }

    i64 powerful = 0;
    i64 squarefree_count = 0;
    for (u64 b = 1; b <= vmax; ++b) {
        if (!sqfree[(size_t)b]) continue;
        ++squarefree_count;
        const u128 b3 = (u128)b * b * b;
        const u64 q = (u64)((u128)N / b3);
        powerful += (i64)isqrt_u64(q);
    }

    const u64 sqrtN = isqrt_u64(N);
    const i64 cubefree_count = count_cubefree_upto(sqrtN);

    const int r6 = (int)iroot6_u64(N);
    const i64 prime6 = (i64)primes_upto(r6).size();

    return powerful - squarefree_count - cubefree_count - prime6 + 1;
}

int main() {
    assert(F(100) == 2);
    assert(F(20'000) == 130);
    assert(F(3'000'000) == 2014);
    std::cout << F(9'000'000'000'000'000'000ULL) << "\n";
    return 0;
}

