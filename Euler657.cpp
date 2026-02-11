#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % kMod);
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) result = mod_mul(result, base);
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 solve_case(int alpha, u64 n) {
    std::vector<u64> inv(static_cast<std::size_t>(alpha + 1), 0);
    inv[1] = 1;
    for (int i = 2; i <= alpha; ++i) {
        inv[static_cast<std::size_t>(i)] =
            kMod - mod_mul(kMod / static_cast<u64>(i), inv[static_cast<std::size_t>(kMod % static_cast<u64>(i))]);
    }

    u64 ans = 0;
    u64 comb = 1;
    const u64 n1 = (n + 1ULL) % kMod;
    const u64 exp = n + 1ULL;

    for (int k = 0; k <= alpha - 1; ++k) {
        u64 geo = 0;
        if (k == 0) {
            geo = 1;
        } else if (k == 1) {
            geo = n1;
        } else {
            const u64 p = mod_pow(static_cast<u64>(k), exp);
            const u64 num = (p + kMod - 1) % kMod;
            geo = mod_mul(num, inv[static_cast<std::size_t>(k - 1)]);
        }

        const u64 term = mod_mul(comb, geo);
        const bool positive = (((alpha - k + 1) & 1) == 0);
        if (positive) {
            ans += term;
            if (ans >= kMod) ans -= kMod;
        } else {
            ans = (ans + kMod - term) % kMod;
        }

        if (k < alpha - 1) {
            comb = mod_mul(comb, static_cast<u64>(alpha - k));
            comb = mod_mul(comb, inv[static_cast<std::size_t>(k + 1)]);
        }
    }

    return ans;
}

}  // namespace

int main() {
    assert(solve_case(3, 0) == 1);
    assert(solve_case(3, 2) == 13);
    assert(solve_case(3, 4) == 79);

    std::cout << solve_case(10'000'000, 1'000'000'000'000ULL) << "\n";
    return 0;
}
