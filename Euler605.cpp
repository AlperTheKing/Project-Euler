#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>

// Project Euler 605: waiting for the first "RL" in a fair coin sequence.
// P_n(k) = sum_{j>=0} (k+jn-1)/2^{k+jn} gives a closed form:
// P_n(k) = 2^{n-k} * ((k-1)2^n + (n-k+1)) / (2^n-1)^2.

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u128 gcd_u128(u128 a, u128 b) {
    while (b) {
        u128 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u128 r = 1;
    u128 x = a % mod;
    while (e) {
        if (e & 1) r = (r * x) % mod;
        x = (x * x) % mod;
        e >>= 1;
    }
    return (u64)r;
}

static u128 M_small_exact(u64 n, u64 k) {
    // Exact only for small n (here used just for statement validations).
    const u128 two_n = (u128)1 << n;
    const u128 two_nk = (u128)1 << (n - k);
    const u128 X = two_n - 1;        // 2^n - 1
    const u128 den = X * X;          // (2^n - 1)^2
    const u128 inner = (u128)(k - 1) * two_n + (u128)(n - k + 1);
    const u128 g = gcd_u128(inner, den);
    const u128 num_red = two_nk * (inner / g);
    const u128 den_red = den / g;
    return num_red * den_red;
}

int main() {
    // Statement validations.
    if ((u64)M_small_exact(3, 1) != 588ULL) {
        std::cerr << "Validation failed: M_3(1)\n";
        return 1;
    }
    if ((u64)M_small_exact(6, 2) != 486864ULL) {
        std::cerr << "Validation failed: M_6(2)\n";
        return 1;
    }

    constexpr u64 MOD = 100000000ULL; // last 8 digits
    constexpr u64 n = 100000007ULL;
    constexpr u64 k = 10007ULL;

    // For this query n is prime, so gcd(n,2^n-1)=1 and the fraction is already reduced.
    const u64 two_n_mod_n = mod_pow(2, n, n);
    const u64 x_mod_n = (two_n_mod_n + n - 1) % n; // (2^n-1) mod n
    if (std::gcd(n, x_mod_n) != 1) {
        std::cerr << "Unexpected: reduction needed (not implemented for large n)\n";
        return 1;
    }

    const u64 two_n_mod = mod_pow(2, n, MOD);
    const u64 two_nk_mod = mod_pow(2, n - k, MOD);
    const u64 x_mod = (two_n_mod + MOD - 1) % MOD;

    const u64 inner_mod = (u64)(((u128)(k - 1) * two_n_mod + (u128)(n - k + 1)) % MOD);

    u64 ans = (u64)(((u128)two_nk_mod * inner_mod) % MOD);
    ans = (u64)(((u128)ans * x_mod) % MOD);
    ans = (u64)(((u128)ans * x_mod) % MOD);

    std::cout << std::setw(8) << std::setfill('0') << ans << "\n";
    return 0;
}

