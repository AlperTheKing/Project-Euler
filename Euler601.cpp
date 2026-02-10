#include <cstdint>
#include <iostream>

// Project Euler 601: Divisibility streaks
//
// For streak(n) >= s we need:
//   n+1 divisible by 2,
//   n+2 divisible by 3,
//   ...
//   n+(s-1) divisible by s,
// which is equivalent to n ≡ 1 (mod m) for all m=2..s.
// Let L_s = lcm(1,2,...,s). Then streak(n) >= s  <=>  n ≡ 1 (mod L_s).
//
// Therefore, for 1 < n < N:
//   P(s,N) = #{n ≡ 1 (mod L_s)} - #{n ≡ 1 (mod L_{s+1})}
//          = floor((N-2)/L_s) - floor((N-2)/L_{s+1}),
// since the solutions are n = 1 + k L with k >= 1 and 1 + kL <= N-1.

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 gcd_u64(u64 a, u64 b) {
    while (b) {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static u128 lcm_u128(u128 a, u64 b) {
    const u64 g = gcd_u64((u64)(a % b), b);
    return a / g * (u128)b;
}

static u128 count_congruent_1(u128 N, u128 L) {
    if (N <= 2 || L == 0) return 0;
    const u128 x = N - 2;
    return (L > x) ? 0 : (x / L);
}

static u128 P(u64 s, u128 N, const u128* L) {
    return count_congruent_1(N, L[s]) - count_congruent_1(N, L[s + 1]);
}

int main() {
    constexpr u64 S_MAX = 32; // need L_{i+1} for i up to 31
    u128 L[S_MAX + 1];
    L[0] = 1;
    for (u64 i = 1; i <= S_MAX; ++i) L[i] = lcm_u128(L[i - 1], i);

    // Statement validations.
    if (P(3, 14, L) != 1) {
        std::cerr << "Validation failed: P(3,14)\n";
        return 1;
    }
    if (P(6, 1000000, L) != 14286) {
        std::cerr << "Validation failed: P(6,10^6)\n";
        return 1;
    }

    u128 ans = 0;
    u128 N = 1;
    for (u64 i = 1; i <= 31; ++i) {
        N *= 4; // N = 4^i
        ans += P(i, N, L);
    }

    // Print u128.
    if (ans == 0) {
        std::cout << "0\n";
        return 0;
    }
    char buf[64];
    int n = 0;
    while (ans > 0) {
        buf[n++] = (char)('0' + (u64)(ans % 10));
        ans /= 10;
    }
    while (n--) std::cout << buf[n];
    std::cout << "\n";
    return 0;
}

