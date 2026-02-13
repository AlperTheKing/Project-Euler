#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    u64 cur = base % kMod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mod_mul(result, cur);
        }
        cur = mod_mul(cur, cur);
        exp >>= 1ULL;
    }
    return result;
}

u64 F_mod(int n) {
    std::vector<u64> fact(static_cast<std::size_t>(n + 1), 1ULL);
    std::vector<u64> inv_fact(static_cast<std::size_t>(n + 1), 1ULL);

    for (int i = 1; i <= n; ++i) {
        fact[static_cast<std::size_t>(i)] = mod_mul(fact[static_cast<std::size_t>(i - 1)], static_cast<u64>(i));
    }
    inv_fact[static_cast<std::size_t>(n)] = mod_pow(fact[static_cast<std::size_t>(n)], kMod - 2ULL);
    for (int i = n; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] = mod_mul(inv_fact[static_cast<std::size_t>(i)], static_cast<u64>(i));
    }

    auto comb = [&](int nn, int kk) -> u64 {
        return mod_mul(fact[static_cast<std::size_t>(nn)],
                       mod_mul(inv_fact[static_cast<std::size_t>(kk)], inv_fact[static_cast<std::size_t>(nn - kk)]));
    };

    u64 sum = 0ULL;
    for (int k = 1; k < n; ++k) {
        const u64 w = static_cast<u64>(std::min(k, n - k));
        u64 term = comb(n, k);
        term = mod_mul(term, w);
        term = mod_mul(term, mod_pow(static_cast<u64>(k), static_cast<u64>(k - 1)));
        term = mod_mul(term, mod_pow(static_cast<u64>(n - k), static_cast<u64>(n - k - 1)));
        sum += term;
        if (sum >= kMod) {
            sum -= kMod;
        }
    }

    const u64 inv2 = (kMod + 1ULL) / 2ULL;
    sum = mod_mul(sum, inv2);
    return mod_mul(fact[static_cast<std::size_t>(n - 1)], sum);
}

u64 pow_u64(u64 base, int exp) {
    u64 result = 1ULL;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

u64 comb_u64(int n, int k) {
    if (k < 0 || k > n) {
        return 0ULL;
    }
    if (k > n - k) {
        k = n - k;
    }
    u64 num = 1ULL;
    u64 den = 1ULL;
    for (int i = 1; i <= k; ++i) {
        num *= static_cast<u64>(n - k + i);
        den *= static_cast<u64>(i);
    }
    return num / den;
}

u64 F_exact_small(int n) {
    u128 sum2 = 0;
    for (int k = 1; k < n; ++k) {
        const u64 w = static_cast<u64>(std::min(k, n - k));
        const u64 term =
            comb_u64(n, k) * w * pow_u64(static_cast<u64>(k), k - 1) * pow_u64(static_cast<u64>(n - k), n - k - 1);
        sum2 += term;
    }
    const u128 sum = sum2 / 2;

    u128 fact = 1;
    for (int i = 2; i <= n - 1; ++i) {
        fact *= static_cast<u64>(i);
    }
    return static_cast<u64>(fact * sum);
}

void run_validations() {
    assert(F_exact_small(3) == 12ULL);
    assert(F_exact_small(4) == 360ULL);
    assert(F_exact_small(8) == 16'785'941'760ULL);
}

}  // namespace

int main() {
    run_validations();
    std::cout << F_mod(100) << '\n';
    return 0;
}
