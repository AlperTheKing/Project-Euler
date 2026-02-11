#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * base) % kMod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % kMod);
        exp >>= 1ULL;
    }
    return result;
}

void batch_invert_range(const u64 left,
                        const u64 right,
                        std::vector<u64>& inverses) {
    const std::size_t m = static_cast<std::size_t>(right - left + 1ULL);
    inverses.assign(m, 0ULL);

    std::vector<u64> prefix(m, 0ULL);
    for (std::size_t i = 0; i < m; ++i) {
        const u64 value = (left + static_cast<u64>(i)) % kMod;
        prefix[i] = (i == 0U)
                        ? value
                        : static_cast<u64>((static_cast<u128>(prefix[i - 1U]) * value) % kMod);
    }

    u64 suffix_inv = mod_pow(prefix[m - 1U], kMod - 2ULL);
    for (std::size_t i = m; i-- > 0U;) {
        const u64 left_prod = (i == 0U) ? 1ULL : prefix[i - 1U];
        inverses[i] = static_cast<u64>((static_cast<u128>(suffix_inv) * left_prod) % kMod);

        const u64 value = (left + static_cast<u64>(i)) % kMod;
        suffix_inv = static_cast<u64>((static_cast<u128>(suffix_inv) * value) % kMod);
    }
}

u64 A(const u64 k, const u64 n) {
    assert(n % k == 0ULL);

    const u64 q = n / k;
    const u64 b = mod_pow(2ULL, q);

    u64 term = mod_pow(b, k);
    u64 answer = term;

    const u64 inv_b = mod_pow(b, kMod - 2ULL);
    const u64 inv_b2 = static_cast<u64>((static_cast<u128>(inv_b) * inv_b) % kMod);

    const u64 half = k / 2ULL;
    constexpr u64 kBlock = 1'000'000ULL;

    std::vector<u64> inverses;

    u64 processed = 0ULL;
    while (processed < half) {
        const u64 left = processed + 1ULL;
        const u64 right = std::min(half, processed + kBlock);

        batch_invert_range(left, right, inverses);

        for (u64 x = left; x <= right; ++x) {
            const u64 inv_x = inverses[static_cast<std::size_t>(x - left)];
            const u64 num1 = (k - 2ULL * (x - 1ULL)) % kMod;
            const u64 num2 = (k - 2ULL * (x - 1ULL) - 1ULL) % kMod;

            term = static_cast<u64>((static_cast<u128>(term) * num1) % kMod);
            term = static_cast<u64>((static_cast<u128>(term) * num2) % kMod);
            term = static_cast<u64>((static_cast<u128>(term) * inv_x) % kMod);
            term = static_cast<u64>((static_cast<u128>(term) * inv_x) % kMod);
            term = static_cast<u64>((static_cast<u128>(term) * inv_b2) % kMod);

            answer += term;
            if (answer >= kMod) {
                answer -= kMod;
            }
        }

        processed = right;
    }

    return answer;
}

}  // namespace

int main() {
    assert(A(3ULL, 9ULL) == 560ULL);
    assert(A(4ULL, 20ULL) == 1'060'870ULL);

    std::cout << A(100'000'000ULL, 10'000'000'000'000'000ULL) << '\n';
    return 0;
}
