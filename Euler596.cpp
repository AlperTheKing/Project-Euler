#include <cmath>
#include <cstdint>
#include <iostream>

// Project Euler 596: Number of Lattice Points in a Hyperball
//
// Let N = r^2. Then
//   T(r) = sum_{n=0}^{N} r_4(n),
// where r_4(n) is the number of representations of n as a sum of four squares.
// For n>0, Jacobi's four-square theorem gives:
//   r_4(n) = 8 * sum_{d|n, 4 \nmid d} d.
// With r_4(0)=1, we obtain
//   T(r) = 1 + 8 * sum_{n=1}^{N} sum_{d|n, 4\nmid d} d
//        = 1 + 8 * sum_{d<=N, 4\nmid d} d * floor(N/d).
// Let
//   S(N) = sum_{d=1}^{N} d * floor(N/d).
// Removing d divisible by 4 yields
//   sum_{d<=N,4\nmid d} d*floor(N/d) = S(N) - sum_{k<=N/4} (4k)*floor(N/(4k))
//                                  = S(N) - 4*S(N/4).
// Therefore
//   T(r) = 1 + 8 * (S(N) - 4*S(N/4)).
//
// We need T(1e8) mod 1e9+7 with N=1e16. We compute S(N) mod M in O(sqrt N) using a
// split hyperbola method: handle i<=sqrt(N) directly, and for i>sqrt(N) group by q=floor(N/i).

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 MOD = 1000000007ULL;
static constexpr u64 INV2 = 500000004ULL; // (MOD+1)/2

static inline u64 mod_add(u64 a, u64 b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}

static inline u64 mod_sub(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + MOD - b);
}

static inline u64 mod_mul(u64 a, u64 b) {
    return (u64)((u128)a * (u128)b % (u128)MOD);
}

static u64 isqrt_u64(u64 x) {
    u64 r = (u64)std::sqrt((long double)x);
    while ((u128)(r + 1) * (u128)(r + 1) <= (u128)x) ++r;
    while ((u128)r * (u128)r > (u128)x) --r;
    return r;
}

static inline u64 sum_arith_mod(u64 l, u64 r) {
    // sum_{i=l}^r i mod MOD
    const u64 cnt = (r - l + 1) % MOD;
    const u64 lr = (l % MOD + r % MOD) % MOD;
    return mod_mul(mod_mul(lr, cnt), INV2);
}

static u64 S_mod(u64 N) {
    if (N == 0) return 0;

    const u64 M = isqrt_u64(N);
    u64 ans = 0;

    // i <= M
    for (u64 i = 1; i <= M; ++i) {
        const u64 q = N / i;
        ans = mod_add(ans, mod_mul(i % MOD, q % MOD));
    }

    // i > M, grouped by q = floor(N/i)
    const u64 qmax = N / (M + 1);
    for (u64 q = 1; q <= qmax; ++q) {
        u64 l = N / (q + 1) + 1;
        const u64 r = N / q;
        if (r <= M) continue;
        if (l <= M) l = M + 1;
        if (l > r) continue;
        ans = mod_add(ans, mod_mul(q % MOD, sum_arith_mod(l, r)));
    }

    return ans;
}

static u64 T_mod(u64 r) {
    const u64 N = r * r;
    const u64 sN = S_mod(N);
    const u64 sN4 = S_mod(N / 4);

    const u64 term = mod_sub(sN, mod_mul(4 % MOD, sN4));
    const u64 t = mod_add(1, mod_mul(8 % MOD, term));
    return t;
}

int main() {
    // Validation points from the statement.
    if (T_mod(2) != 89ULL) {
        std::cerr << "Validation failed: T(2)\n";
        return 1;
    }
    if (T_mod(5) != 3121ULL) {
        std::cerr << "Validation failed: T(5)\n";
        return 1;
    }
    if (T_mod(100) != 493490641ULL) {
        std::cerr << "Validation failed: T(100)\n";
        return 1;
    }
    {
        const u64 expected = 49348022079085897ULL % MOD;
        if (T_mod(10000) != expected) {
            std::cerr << "Validation failed: T(1e4)\n";
            return 1;
        }
    }

    std::cout << T_mod(100000000ULL) << "\n";
    return 0;
}
