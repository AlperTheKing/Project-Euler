#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'234'567'891ULL;
constexpr u64 kSplit = 3'000'000ULL;

inline u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= kMod) {
        a -= kMod;
    }
    return a;
}

inline u64 sub_mod(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

inline u64 mul_mod(u64 a, u64 b) {
    return (a * b) % kMod;
}

u64 pow_mod(u64 base, u64 exp) {
    u64 result = 1 % kMod;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 sum_range_mod(u64 l, u64 r) {
    if (l > r) {
        return 0;
    }
    const u128 cnt = static_cast<u128>(r) - static_cast<u128>(l) + 1U;
    const u128 s = static_cast<u128>(l) + static_cast<u128>(r);
    const u128 total = (cnt * s) / 2U;
    return static_cast<u64>(total % kMod);
}

u64 solve(u64 n) {
    const u64 pow_n_minus_1 = pow_mod(2, n - 1);
    const u64 split = (n < kSplit) ? n : kSplit;

    u64 weighted_sum = 0;

    for (u64 x = 1; x <= split; ++x) {
        const u64 q = n / x;
        const u64 p = pow_mod(2, n - q);
        weighted_sum = add_mod(weighted_sum, mul_mod(x % kMod, p));
    }

    const u64 inv2 = (kMod + 1ULL) / 2ULL;
    const u64 qmax = n / (split + 1ULL);
    u64 p = pow_n_minus_1;  // 2^(n-1)

    for (u64 q = 1; q <= qmax; ++q) {
        u64 l = n / (q + 1ULL) + 1ULL;
        if (l <= split) {
            l = split + 1ULL;
        }
        const u64 r = n / q;
        if (l <= r) {
            const u64 sum_x = sum_range_mod(l, r);
            weighted_sum = add_mod(weighted_sum, mul_mod(sum_x, p));
        }
        p = mul_mod(p, inv2);  // 2^(n-(q+1))
    }

    const u64 tri = sum_range_mod(1, n);
    const u64 first = mul_mod(pow_n_minus_1, tri);
    return sub_mod(first, weighted_sum);
}

u64 solve_naive(u64 n) {
    u64 ans = 0;
    const u64 pow_n_minus_1 = pow_mod(2, n - 1);
    for (u64 x = 1; x <= n; ++x) {
        const u64 q = n / x;
        const u64 ways = sub_mod(pow_n_minus_1, pow_mod(2, n - q));
        ans = add_mod(ans, mul_mod(x % kMod, ways));
    }
    return ans;
}

void run_validations() {
    assert(solve(1) == 0ULL);
    assert(solve(10) == 4'927ULL);
    for (u64 n = 2; n <= 200; ++n) {
        assert(solve(n) == solve_naive(n));
    }
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kN = 100'000'000'000'000ULL;
    std::cout << solve(kN) << '\n';
    return 0;
}
