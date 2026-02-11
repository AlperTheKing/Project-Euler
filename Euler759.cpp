#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;
using u64 = std::uint64_t;

constexpr i64 MOD = 1'000'000'007LL;
constexpr int MAX_BITS = 64;

struct Aggregates {
    std::array<i64, MAX_BITS + 1> cnt{};
    std::array<i64, MAX_BITS + 1> sx{};
    std::array<i64, MAX_BITS + 1> sx2{};
    std::array<i64, MAX_BITS + 1> sc{};
    std::array<i64, MAX_BITS + 1> sc2{};
    std::array<i64, MAX_BITS + 1> sxc{};
    std::array<i64, MAX_BITS + 1> sxc2{};
    std::array<i64, MAX_BITS + 1> sx2c{};
    std::array<i64, MAX_BITS + 1> sx2c2{};
};

i64 add_mod(i64 a, i64 b) {
    i64 r = a + b;
    if (r >= MOD) {
        r -= MOD;
    }
    return r;
}

i64 sub_mod(i64 a, i64 b) {
    i64 r = a - b;
    if (r < 0) {
        r += MOD;
    }
    return r;
}

i64 mul_mod(i64 a, i64 b) {
    return static_cast<i64>((static_cast<i128>(a) * b) % MOD);
}

Aggregates build_aggregates() {
    Aggregates ag;
    ag.cnt[0] = 1;

    for (int m = 0; m < MAX_BITS; ++m) {
        const i64 n = ag.cnt[m];
        const i64 sx = ag.sx[m];
        const i64 sx2 = ag.sx2[m];
        const i64 sc = ag.sc[m];
        const i64 sc2 = ag.sc2[m];
        const i64 sxc = ag.sxc[m];
        const i64 sxc2 = ag.sxc2[m];
        const i64 sx2c = ag.sx2c[m];
        const i64 sx2c2 = ag.sx2c2[m];

        const i64 two_sx = mul_mod(2, sx);
        const i64 four_sx = mul_mod(4, sx);
        const i64 four_sx2 = mul_mod(4, sx2);
        const i64 two_sc = mul_mod(2, sc);

        const i64 e_cnt = n;
        const i64 e_sx = two_sx;
        const i64 e_sx2 = four_sx2;
        const i64 e_sc = sc;
        const i64 e_sc2 = sc2;
        const i64 e_sxc = mul_mod(2, sxc);
        const i64 e_sxc2 = mul_mod(2, sxc2);
        const i64 e_sx2c = mul_mod(4, sx2c);
        const i64 e_sx2c2 = mul_mod(4, sx2c2);

        const i64 o_cnt = n;
        const i64 o_sx = add_mod(two_sx, n);
        const i64 o_sx2 = add_mod(add_mod(four_sx2, four_sx), n);
        const i64 o_sc = add_mod(sc, n);
        const i64 o_sc2 = add_mod(add_mod(sc2, two_sc), n);

        i64 o_sxc = 0;
        o_sxc = add_mod(o_sxc, mul_mod(2, sxc));
        o_sxc = add_mod(o_sxc, sc);
        o_sxc = add_mod(o_sxc, two_sx);
        o_sxc = add_mod(o_sxc, n);

        i64 o_sxc2 = 0;
        o_sxc2 = add_mod(o_sxc2, mul_mod(2, sxc2));
        o_sxc2 = add_mod(o_sxc2, sc2);
        o_sxc2 = add_mod(o_sxc2, mul_mod(4, sxc));
        o_sxc2 = add_mod(o_sxc2, mul_mod(2, sc));
        o_sxc2 = add_mod(o_sxc2, two_sx);
        o_sxc2 = add_mod(o_sxc2, n);

        i64 o_sx2c = 0;
        o_sx2c = add_mod(o_sx2c, mul_mod(4, sx2c));
        o_sx2c = add_mod(o_sx2c, mul_mod(4, sxc));
        o_sx2c = add_mod(o_sx2c, sc);
        o_sx2c = add_mod(o_sx2c, four_sx2);
        o_sx2c = add_mod(o_sx2c, four_sx);
        o_sx2c = add_mod(o_sx2c, n);

        i64 o_sx2c2 = 0;
        o_sx2c2 = add_mod(o_sx2c2, mul_mod(4, sx2c2));
        o_sx2c2 = add_mod(o_sx2c2, mul_mod(4, sxc2));
        o_sx2c2 = add_mod(o_sx2c2, sc2);
        o_sx2c2 = add_mod(o_sx2c2, mul_mod(8, sx2c));
        o_sx2c2 = add_mod(o_sx2c2, mul_mod(8, sxc));
        o_sx2c2 = add_mod(o_sx2c2, mul_mod(2, sc));
        o_sx2c2 = add_mod(o_sx2c2, four_sx2);
        o_sx2c2 = add_mod(o_sx2c2, four_sx);
        o_sx2c2 = add_mod(o_sx2c2, n);

        ag.cnt[m + 1] = add_mod(e_cnt, o_cnt);
        ag.sx[m + 1] = add_mod(e_sx, o_sx);
        ag.sx2[m + 1] = add_mod(e_sx2, o_sx2);
        ag.sc[m + 1] = add_mod(e_sc, o_sc);
        ag.sc2[m + 1] = add_mod(e_sc2, o_sc2);
        ag.sxc[m + 1] = add_mod(e_sxc, o_sxc);
        ag.sxc2[m + 1] = add_mod(e_sxc2, o_sxc2);
        ag.sx2c[m + 1] = add_mod(e_sx2c, o_sx2c);
        ag.sx2c2[m + 1] = add_mod(e_sx2c2, o_sx2c2);
    }

    return ag;
}

i64 block_sum(i64 prefix_value_mod, i64 prefix_popcount, int lower_bits, const Aggregates& ag) {
    const i64 p = prefix_value_mod;
    const i64 c = prefix_popcount % MOD;
    const i64 p2 = mul_mod(p, p);
    const i64 c2 = mul_mod(c, c);

    i64 ans = 0;
    ans = add_mod(ans, mul_mod(mul_mod(ag.cnt[lower_bits], p2), c2));
    ans = add_mod(ans, mul_mod(ag.sc[lower_bits], mul_mod(mul_mod(2, p2), c)));
    ans = add_mod(ans, mul_mod(ag.sc2[lower_bits], p2));
    ans = add_mod(ans, mul_mod(ag.sx[lower_bits], mul_mod(mul_mod(2, p), c2)));
    ans = add_mod(ans, mul_mod(ag.sxc[lower_bits], mul_mod(mul_mod(4, p), c)));
    ans = add_mod(ans, mul_mod(ag.sxc2[lower_bits], mul_mod(2, p)));
    ans = add_mod(ans, mul_mod(ag.sx2[lower_bits], c2));
    ans = add_mod(ans, mul_mod(ag.sx2c[lower_bits], mul_mod(2, c)));
    ans = add_mod(ans, ag.sx2c2[lower_bits]);
    return ans;
}

i64 S(u64 n, const Aggregates& ag, const std::array<i64, MAX_BITS + 1>& pow2_mod) {
    i64 ans = 0;
    i64 prefix_value_mod = 0;
    i64 prefix_popcount = 0;

    for (int bit = MAX_BITS - 1; bit >= 0; --bit) {
        if (((n >> bit) & 1ULL) == 0ULL) {
            continue;
        }
        ans = add_mod(ans, block_sum(prefix_value_mod, prefix_popcount, bit, ag));
        prefix_value_mod = add_mod(prefix_value_mod, pow2_mod[bit]);
        ++prefix_popcount;
    }

    const i64 p2 = mul_mod(prefix_value_mod, prefix_value_mod);
    const i64 c = prefix_popcount % MOD;
    const i64 c2 = mul_mod(c, c);
    ans = add_mod(ans, mul_mod(p2, c2));
    return ans;
}

u64 brute_S(int n) {
    std::vector<u64> f(static_cast<std::size_t>(n + 1), 0ULL);
    f[1] = 1ULL;
    for (int i = 2; i <= n; ++i) {
        if ((i & 1) == 0) {
            f[i] = 2ULL * f[i / 2];
        } else {
            int m = (i - 1) / 2;
            f[i] = static_cast<u64>(2 * m + 1) + 2ULL * f[m] + f[m] / static_cast<u64>(m);
        }
    }
    u64 sum = 0ULL;
    for (int i = 1; i <= n; ++i) {
        sum += f[i] * f[i];
    }
    return sum;
}

}  // namespace

int main() {
    const Aggregates ag = build_aggregates();

    std::array<i64, MAX_BITS + 1> pow2_mod{};
    pow2_mod[0] = 1;
    for (int i = 1; i <= MAX_BITS; ++i) {
        pow2_mod[i] = mul_mod(2, pow2_mod[i - 1]);
    }

    assert(brute_S(10) == 1530ULL);
    assert(brute_S(100) == 4'798'445ULL);
    assert(S(10ULL, ag, pow2_mod) == 1530LL);
    assert(S(100ULL, ag, pow2_mod) == 4'798'445LL);

    std::cout << S(10'000'000'000'000'000ULL, ag, pow2_mod) << '\n';
    return 0;
}
