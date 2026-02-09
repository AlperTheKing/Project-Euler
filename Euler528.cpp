#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_pow(u64 a, u64 e) {
    u64 r = 1 % kMod;
    a %= kMod;
    while (e) {
        if (e & 1) {
            r = static_cast<u128>(r) * a % kMod;
        }
        a = static_cast<u128>(a) * a % kMod;
        e >>= 1;
    }
    return r;
}

u64 binom_large_n_small_k(const u128 n, const int k, const std::vector<u64>& inv_fact) {
    if (k < 0) {
        return 0;
    }
    const u64 nm = static_cast<u64>(n % static_cast<u128>(kMod));
    if (nm < static_cast<u64>(k)) {
        return 0;
    }
    u64 num = 1;
    for (int i = 0; i < k; ++i) {
        num = static_cast<u128>(num) * (nm - static_cast<u64>(i)) % kMod;
    }
    return static_cast<u128>(num) * inv_fact[static_cast<std::size_t>(k)] % kMod;
}

u64 S(const u128 n, const int k, const int b, const std::vector<u64>& inv_fact) {
    // Count solutions to x1+...+xk <= n with 0<=xm<=b^m using inclusion-exclusion over an added slack variable:
    // x1+...+xk + s = n, s>=0, and upper bounds via shifts by (b^m+1).
    const int masks = 1 << k;

    std::vector<u128> add(static_cast<std::size_t>(k));
    u128 p = 1;
    for (int m = 1; m <= k; ++m) {
        p *= static_cast<u128>(b);
        add[static_cast<std::size_t>(m - 1)] = p + 1;  // (b^m + 1)
    }

    std::vector<u128> shift(static_cast<std::size_t>(masks), 0);
    for (int mask = 1; mask < masks; ++mask) {
        const int bit = __builtin_ctz(static_cast<unsigned>(mask));
        const int prev = mask & (mask - 1);
        shift[static_cast<std::size_t>(mask)] =
            shift[static_cast<std::size_t>(prev)] + add[static_cast<std::size_t>(bit)];
    }

    u64 ans = 0;
    for (int mask = 0; mask < masks; ++mask) {
        const u128 sh = shift[static_cast<std::size_t>(mask)];
        if (sh > n) {
            continue;
        }
        const u128 N = (n - sh) + static_cast<u128>(k);
        const u64 term = binom_large_n_small_k(N, k, inv_fact);
        if ((__builtin_popcount(static_cast<unsigned>(mask)) & 1) == 0) {
            ans += term;
            if (ans >= kMod) {
                ans -= kMod;
            }
        } else {
            ans = (ans >= term) ? (ans - term) : (ans + kMod - term);
        }
    }
    return ans;
}

u128 pow10_u128(int e) {
    u128 r = 1;
    for (int i = 0; i < e; ++i) {
        r *= 10;
    }
    return r;
}

bool run_checkpoints() {
    // inv_fact up to 15 is sufficient for all tasks here.
    std::vector<u64> fact(16, 1), inv_fact(16, 1);
    for (int i = 1; i <= 15; ++i) {
        fact[static_cast<std::size_t>(i)] = static_cast<u128>(fact[static_cast<std::size_t>(i - 1)]) * i % kMod;
    }
    inv_fact[15] = mod_pow(fact[15], kMod - 2);
    for (int i = 15; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<u128>(inv_fact[static_cast<std::size_t>(i)]) * i % kMod;
    }

    if (S(14, 3, 2, inv_fact) != 135ULL) {
        std::cerr << "Checkpoint failed: S(14,3,2)\n";
        return false;
    }
    if (S(200, 5, 3, inv_fact) != 12'949'440ULL) {
        std::cerr << "Checkpoint failed: S(200,5,3)\n";
        return false;
    }
    if (S(1000, 10, 5, inv_fact) != 624'839'075ULL) {
        std::cerr << "Checkpoint failed: S(1000,10,5) mod\n";
        return false;
    }
    return true;
}

u64 solve() {
    std::vector<u64> fact(16, 1), inv_fact(16, 1);
    for (int i = 1; i <= 15; ++i) {
        fact[static_cast<std::size_t>(i)] = static_cast<u128>(fact[static_cast<std::size_t>(i - 1)]) * i % kMod;
    }
    inv_fact[15] = mod_pow(fact[15], kMod - 2);
    for (int i = 15; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<u128>(inv_fact[static_cast<std::size_t>(i)]) * i % kMod;
    }

    u64 total = 0;
    for (int k = 10; k <= 15; ++k) {
        const u64 term = S(pow10_u128(k), k, k, inv_fact);
        total += term;
        total %= kMod;
    }
    return total;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }
    std::cout << solve() << '\n';
    return 0;
}

