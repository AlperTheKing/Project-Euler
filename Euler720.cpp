#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 S_pow2_mod(const int exponent) {
    if (exponent == 1) {
        return 1;
    }
    if (exponent == 2) {
        return 3;
    }

    std::vector<u32> perm{1, 3, 2, 4};
    std::vector<u32> lehmer{0, 1, 0, 0};

    u32 n = 4;
    const u32 target = 1u << exponent;

    while (n < target) {
        perm.resize(static_cast<std::size_t>(2 * n));
        lehmer.resize(static_cast<std::size_t>(2 * n));

        for (int i = static_cast<int>(n) - 1; i >= 1; --i) {
            const std::size_t idx = static_cast<std::size_t>(i);
            const std::size_t out = static_cast<std::size_t>(n) + idx;
            perm[out] = 2u * perm[idx];
            lehmer[out] = lehmer[idx];
        }

        perm[static_cast<std::size_t>(n)] = 2u * perm[static_cast<std::size_t>(n - 1)] - 1u;
        lehmer[static_cast<std::size_t>(n)] = n - 2u;

        perm[static_cast<std::size_t>(n - 1)] = 2u * perm[0];
        lehmer[static_cast<std::size_t>(n - 1)] = 0u;

        for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
            const std::size_t idx = static_cast<std::size_t>(i);
            lehmer[idx] = perm[idx] - 1u + lehmer[idx];
            perm[idx] = 2u * perm[idx] - 1u;
        }

        n *= 2u;
    }

    u64 rank = 1;
    u64 fact = 1;
    for (int i = static_cast<int>(target) - 1, k = 0; i >= 0; --i, ++k) {
        rank = (rank + static_cast<u64>(lehmer[static_cast<std::size_t>(i)]) * fact) % kMod;
        fact = (fact * static_cast<u64>(k + 1)) % kMod;
    }

    return rank;
}

}  // namespace

int main() {
    assert(S_pow2_mod(2) == 3);
    assert(S_pow2_mod(3) == 2295);
    assert(S_pow2_mod(5) == 641'839'205ULL);

    std::cout << S_pow2_mod(25) << '\n';
    return 0;
}
