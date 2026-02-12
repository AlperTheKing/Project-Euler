#include <cassert>
#include <cstdint>
#include <iostream>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static inline u64 pow_mod10(u64 exp, u64 mod) {
    u64 base = 10 % mod;
    u64 res = 1 % mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = static_cast<u64>((static_cast<u128>(res) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return res;
}

static inline int nth_digit_reciprocal(u64 n, u64 k) {
    const u64 r = pow_mod10(n - 1, k);
    return static_cast<int>((10ULL * r) / k);
}

static u64 S(u64 n) {
    u64 ans = 0;
    for (u64 k = 1; k <= n; ++k) {
        ans += static_cast<u64>(nth_digit_reciprocal(n, k));
    }
    return ans;
}

int main() {
    assert(S(7) == 10);
    assert(S(100) == 418);
    std::cout << S(10'000'000) << '\n';
    return 0;
}
