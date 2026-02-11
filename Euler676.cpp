#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 MOD16 = 10'000'000'000'000'000ULL;

u64 add_mod(u64 a, u64 b) {
    const u64 s = a + b;
    if (s >= MOD16 || s < a) {
        return s - MOD16;
    }
    return s;
}

u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % MOD16);
}

u64 M(u64 n, int k, int l) {
    if (n == 0ULL) {
        return 0ULL;
    }

    const int max_bit = 64 - __builtin_clzll(n);

    std::vector<int> coeff(static_cast<std::size_t>(max_bit), 0);
    int max_abs = 0;
    for (int p = 0; p < max_bit; ++p) {
        coeff[static_cast<std::size_t>(p)] =
            (1 << (p % k)) - (1 << (p % l));
        max_abs += std::abs(coeff[static_cast<std::size_t>(p)]);
    }

    const int offset = max_abs;
    const int width = 2 * max_abs + 1;

    std::vector<u64> cnt_lo(static_cast<std::size_t>(width), 0ULL);
    std::vector<u64> cnt_hi(static_cast<std::size_t>(width), 0ULL);
    std::vector<u64> sum_lo(static_cast<std::size_t>(width), 0ULL);
    std::vector<u64> sum_hi(static_cast<std::size_t>(width), 0ULL);

    cnt_hi[static_cast<std::size_t>(offset)] = 1ULL;

    for (int pos = max_bit - 1; pos >= 0; --pos) {
        std::vector<u64> ncnt_lo(static_cast<std::size_t>(width), 0ULL);
        std::vector<u64> ncnt_hi(static_cast<std::size_t>(width), 0ULL);
        std::vector<u64> nsum_lo(static_cast<std::size_t>(width), 0ULL);
        std::vector<u64> nsum_hi(static_cast<std::size_t>(width), 0ULL);

        const int lim = static_cast<int>((n >> pos) & 1ULL);
        const int c = coeff[static_cast<std::size_t>(pos)];
        const u64 bit_value = (1ULL << pos) % MOD16;

        for (int idx = 0; idx < width; ++idx) {
            const u64 c0 = cnt_lo[static_cast<std::size_t>(idx)];
            if (c0 != 0ULL) {
                const u64 s0 = sum_lo[static_cast<std::size_t>(idx)];

                ncnt_lo[static_cast<std::size_t>(idx)] += c0;
                nsum_lo[static_cast<std::size_t>(idx)] =
                    add_mod(nsum_lo[static_cast<std::size_t>(idx)], s0);

                const int j = idx + c;
                ncnt_lo[static_cast<std::size_t>(j)] += c0;
                const u64 add1 = add_mod(s0, mul_mod(c0, bit_value));
                nsum_lo[static_cast<std::size_t>(j)] =
                    add_mod(nsum_lo[static_cast<std::size_t>(j)], add1);
            }

            const u64 c1 = cnt_hi[static_cast<std::size_t>(idx)];
            if (c1 == 0ULL) {
                continue;
            }
            const u64 s1 = sum_hi[static_cast<std::size_t>(idx)];

            if (lim == 0) {
                ncnt_hi[static_cast<std::size_t>(idx)] += c1;
                nsum_hi[static_cast<std::size_t>(idx)] =
                    add_mod(nsum_hi[static_cast<std::size_t>(idx)], s1);
            } else {
                ncnt_lo[static_cast<std::size_t>(idx)] += c1;
                nsum_lo[static_cast<std::size_t>(idx)] =
                    add_mod(nsum_lo[static_cast<std::size_t>(idx)], s1);

                const int j = idx + c;
                ncnt_hi[static_cast<std::size_t>(j)] += c1;
                const u64 add1 = add_mod(s1, mul_mod(c1, bit_value));
                nsum_hi[static_cast<std::size_t>(j)] =
                    add_mod(nsum_hi[static_cast<std::size_t>(j)], add1);
            }
        }

        cnt_lo.swap(ncnt_lo);
        cnt_hi.swap(ncnt_hi);
        sum_lo.swap(nsum_lo);
        sum_hi.swap(nsum_hi);
    }

    return add_mod(sum_lo[static_cast<std::size_t>(offset)],
                   sum_hi[static_cast<std::size_t>(offset)]);
}

}  // namespace

int main() {
    assert(M(10ULL, 3, 1) == 18ULL);
    assert(M(100ULL, 3, 1) == 292ULL);
    assert(M(1'000'000ULL, 3, 1) == 19'173'952ULL);

    u64 ans = 0ULL;
    for (int k = 3; k <= 6; ++k) {
        for (int l = 1; l <= k - 2; ++l) {
            ans = add_mod(ans, M(10'000'000'000'000'000ULL, k, l));
        }
    }

    std::cout << std::setw(16) << std::setfill('0') << ans << "\n";
    return 0;
}
