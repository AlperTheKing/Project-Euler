#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u64 v2_u64(u64 x) {
    if (x == 0) return 64;
    return static_cast<u64>(__builtin_ctzll(x));
}

// For p=2, Legendre gives v2(n!) = sum_{k>=1} floor(n/2^k) = n - popcount(n).
static u64 v2_factorial(u64 n) {
    return n - static_cast<u64>(__builtin_popcountll(n));
}

// Directly compute S((p_m#)^n) for small m,n using the divisor-structure formula:
// S(N) = sum_{d|N} tau(d) - tau(N), and for N = (p_m#)^n:
//   tau(N) = (n+1)^m
//   sum_{d|N} tau(d) = (sum_{k=0..n} (k+1))^m = ((n+1)(n+2)/2)^m
static u64 S_primorial_power_small(u64 m, u64 n) {
    const u64 t = (n + 1) * (n + 2) / 2;
    u128 A = 1, B = 1;
    for (u64 i = 0; i < m; ++i) {
        A *= static_cast<u128>(t);
        B *= static_cast<u128>(n + 1);
    }
    const u128 S = A - B;
    return static_cast<u64>(S);  // used only for tiny cases where it fits
}

// For odd m, the 2-adic valuation of S((p_m#)^n) collapses to a simple piecewise rule:
// Let m be odd and n>=1.
//   If n is odd:   E(m,n) = m * v2((n+1)/2)
//   If n≡2 (mod4): E(m,n) = 0
//   If 4|n:        E(m,n) = v2(n) - 1
static u64 E_odd_m(u64 m, u64 n) {
    assert((m & 1) == 1);
    if (n & 1) {
        return m * v2_u64((n + 1) / 2);
    }
    if ((n & 3) == 2) return 0;
    return v2_u64(n) - 1;
}

static u128 Q(u64 m, u64 N) {
    assert((m & 1) == 1);

    // Sum E(m,n) for n=1..N using the piecewise rule:
    // - odd n: contribute m*v2((n+1)/2); letting k=(n+1)/2, k=1..floor((N+1)/2) -> m*v2(K!)
    // - n≡2 (mod4): contribute 0
    // - n=4t: contribute v2(4t)-1 = 1+v2(t); sum_{t<=floor(N/4)} (1+v2(t)) = T + v2(T!)
    const u64 K = (N + 1) / 2;
    const u64 T = N / 4;

    const u128 part_odd = static_cast<u128>(m) * static_cast<u128>(v2_factorial(K));
    const u128 part_mult4 = static_cast<u128>(T) + static_cast<u128>(v2_factorial(T));
    return part_odd + part_mult4;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Statement check: E(2,1)=0 since S(6)=5.
    {
        const u64 s = S_primorial_power_small(2, 1);
        assert(s == 5);
        assert(v2_u64(s) == 0);
    }

    // Statement check: Q(8)=2714886 for m=904961.
    {
        const u64 m = 904961;
        u128 sum = 0;
        for (u64 i = 1; i <= 8; ++i) sum += E_odd_m(m, i);
        assert(to_string_u128(sum) == "2714886");
        assert(to_string_u128(Q(m, 8)) == "2714886");
    }

    const u64 m = 904961;
    const u64 N = 1'000'000'000'000ULL;
    std::cout << to_string_u128(Q(m, N)) << '\n';
    return 0;
}

