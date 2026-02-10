#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

// Project Euler 570: Snowflakes
//
// With the correct interpretation of the construction, the counts satisfy:
//   A(n) = 3*4^(n-1) - 2*3^(n-1)
//   B(n) = ((9n-69)*4^(n-1) + 2*(4n+26)*3^(n-1)) / 2
//
// For the gcd, we avoid huge integers. Let
//   x = 4^(n-1), y = 3^(n-1).
// Then
//   A = 3x - 2y
//   2B = (9n-69)x + (8n+52)y
// Any common divisor of A and 2B divides the determinant
//   det = 3*(8n+52) - (-2)*(9n-69) = 42n+18 = 6*(7n+3).
// Empirically and by the derived closed forms, gcd(A,B)=gcd(A,2B)=gcd(A,42n+18).
//
// Also for n>=2:
//   A(n) = 6*(2^(2n-3) - 3^(n-2)) = 6*(2*4^(n-2) - 3^(n-2)).
// Hence
//   G(n) = gcd(A(n), B(n)) = gcd(A(n), 42n+18)
//        = 6 * gcd(7n+3, 2*4^(n-2) - 3^(n-2)).
//
// We compute the latter gcd using modular exponentiation modulo m = 7n+3 (<= 7e7+3),
// which is fast enough for n up to 1e7, and sum G(n).

using u32 = uint32_t;
using u64 = uint64_t;

static inline std::pair<u32, u32> pow4_pow3_mod(u32 exp, u32 mod) {
    // Returns (4^exp mod mod, 3^exp mod mod) using one exponent bit-scan.
    u64 r4 = 1 % mod, r3 = 1 % mod;
    u64 b4 = 4 % mod, b3 = 3 % mod;
    u32 e = exp;
    while (e) {
        if (e & 1U) {
            r4 = (r4 * b4) % mod;
            r3 = (r3 * b3) % mod;
        }
        b4 = (b4 * b4) % mod;
        b3 = (b3 * b3) % mod;
        e >>= 1U;
    }
    return {static_cast<u32>(r4), static_cast<u32>(r3)};
}

static inline u64 G(u32 n) {
    // n >= 3
    const u32 exp = n - 2;
    const u32 m = 7U * n + 3U;
    const auto [p4, p3] = pow4_pow3_mod(exp, m);

    // v = (2*4^exp - 3^exp) mod m, in [0, m).
    u64 v = (2ULL * p4) % m;
    v = (v + m - p3) % m;

    const u32 g = std::gcd(m, static_cast<u32>(v));
    return 6ULL * g;
}

static u64 pow_u64(u64 a, int e) {
    u64 r = 1;
    while (e--) r *= a;
    return r;
}

static u64 A_u64(int n) {
    // Safe for validation ranges in this file (n <= 11).
    return 3ULL * pow_u64(4ULL, n - 1) - 2ULL * pow_u64(3ULL, n - 1);
}

static u64 B_u64(int n) {
    // B(n) = ((9n-69)*4^(n-1) + 2*(4n+26)*3^(n-1))/2
    __int128 x = pow_u64(4ULL, n - 1);
    __int128 y = pow_u64(3ULL, n - 1);
    __int128 num = (__int128(9) * n - 69) * x + __int128(2) * (__int128(4) * n + 26) * y;
    return static_cast<u64>(num / 2);
}

int main() {
    // Validation points from the statement.
    {
        const u64 a3 = A_u64(3), b3 = B_u64(3);
        if (a3 != 30 || b3 != 6 || std::gcd(a3, b3) != 6 || G(3) != 6) {
            std::cerr << "Validation failed at n=3\n";
            return 1;
        }
        const u64 a11 = A_u64(11), b11 = B_u64(11);
        if (a11 != 3027630ULL || b11 != 19862070ULL || std::gcd(a11, b11) != 30 || G(11) != 30) {
            std::cerr << "Validation failed at n=11\n";
            return 1;
        }
        if (G(500) != 186ULL) {
            std::cerr << "Validation failed at n=500\n";
            return 1;
        }
        u64 s = 0;
        for (u32 n = 3; n <= 500; ++n) s += G(n);
        if (s != 5124ULL) {
            std::cerr << "Validation failed for sum_{3..500}\n";
            return 1;
        }
    }

    static constexpr u32 N = 10'000'000;

    const u32 hw = std::thread::hardware_concurrency();
    const u32 threads = std::max<u32>(1U, std::min<u32>(8U, hw ? hw : 4U));

    std::vector<std::thread> pool;
    std::vector<u64> partial(threads, 0);

    const u32 total = N - 3 + 1;
    const u32 block = (total + threads - 1) / threads;

    for (u32 t = 0; t < threads; ++t) {
        const u32 n_lo = 3U + t * block;
        const u32 n_hi = std::min<u32>(N, n_lo + block - 1);
        pool.emplace_back([&, t, n_lo, n_hi] {
            u64 acc = 0;
            for (u32 n = n_lo; n <= n_hi; ++n) acc += G(n);
            partial[t] = acc;
        });
    }
    for (auto& th : pool) th.join();

    u64 ans = 0;
    for (u64 x : partial) ans += x;
    std::cout << ans << "\n";
    return 0;
}

