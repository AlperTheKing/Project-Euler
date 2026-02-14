#include <cassert>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <functional>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr i64 MOD = 1'000'000'007LL;
constexpr i64 INV2 = (MOD + 1) / 2;

u64 isqrt_u64(const u64 n) {
    if (n == 0) {
        return 0;
    }
    long double d = static_cast<long double>(n);
    u64 x = static_cast<u64>(std::sqrt(d));
    while (static_cast<u128>(x + 1) * (x + 1) <= n) {
        ++x;
    }
    while (static_cast<u128>(x) * x > n) {
        --x;
    }
    return x;
}

u64 icbrt_u64(const u64 n) {
    if (n == 0) {
        return 0;
    }
    long double d = static_cast<long double>(n);
    u64 x = static_cast<u64>(std::cbrt(d));
    while (static_cast<u128>(x + 1) * (x + 1) * (x + 1) <= n) {
        ++x;
    }
    while (static_cast<u128>(x) * x * x > n) {
        --x;
    }
    return x;
}

u64 p2_prefix(const u64 r) {
    // Sum_{t=1..r} 2*ceil(2*sqrt(t)); p2_prefix(0)=0.
    if (r == 0) {
        return 0;
    }

    u64 s = isqrt_u64(r);
    if (static_cast<u128>(s) * s < r) {
        ++s;  // s = ceil(sqrt(r))
    }
    const u64 k = s - 1;

    const u128 base = (static_cast<u128>(8) * k * k * k + static_cast<u128>(3) * k * k + k) /
                      static_cast<u128>(3);
    const u64 consumed = k * k;
    const u64 rem = r - consumed;

    const u64 len1 = (rem < (s - 1)) ? rem : (s - 1);
    const u64 len2 = rem - len1;

    const u128 add1 = static_cast<u128>(len1) * (4 * s - 2);
    const u128 add2 = static_cast<u128>(len2) * (4 * s);

    return static_cast<u64>(base + add1 + add2);
}

i64 mod_mul(const i64 a, const i64 b) {
    return static_cast<i64>((static_cast<__int128>(a) * b) % MOD);
}

i64 mod_sum_range(const u64 l, const u64 r) {
    if (l > r) {
        return 0;
    }
    const i64 len = static_cast<i64>((r - l + 1) % MOD);
    const i64 ends = static_cast<i64>((l % MOD + r % MOD) % MOD);
    return mod_mul(mod_mul(len, ends), INV2);
}

i64 add_phase(i64 acc, const u64 L, const u64 R, const u64 capN, const u64 base, const u64 rBase) {
    if (L > capN) {
        return acc;
    }
    const u64 left = L;
    const u64 right = (R < capN) ? R : capN;
    if (left > right) {
        return acc;
    }

    const u64 len_u = right - left + 1;
    const i64 len = static_cast<i64>(len_u % MOD);
    const i64 sum_n = mod_sum_range(left, right);

    const u64 rr = right - rBase;
    const u64 ll = left - rBase;
    const u64 sum_p_u = p2_prefix(rr) - p2_prefix((ll == 0) ? 0 : (ll - 1));
    const i64 sum_p = static_cast<i64>(sum_p_u % MOD);

    i64 cur = mod_mul(6, sum_n);
    cur = (cur - mod_mul(len, static_cast<i64>(base % MOD)) + MOD) % MOD;
    cur = (cur - sum_p + MOD) % MOD;

    acc += cur;
    if (acc >= MOD) {
        acc -= MOD;
    }
    return acc;
}

i64 G_mod(const u64 N) {
    i64 ans = 0;
    const u64 mMax = icbrt_u64(N);

    for (u64 m = 1; m <= mMax; ++m) {
        const u64 m2 = m * m;
        const u64 mp1 = m + 1;

        const u64 v0 = m * m2;
        const u64 v1 = m2 * mp1;
        const u64 v2 = m * mp1 * mp1;
        const u64 v3 = mp1 * mp1 * mp1;

        const u64 s0 = 6 * m2;
        const u64 s1 = 2 * (m2 + 2 * m * mp1);
        const u64 s2 = 2 * (2 * m * mp1 + mp1 * mp1);

        ans = add_phase(ans, v0, v1, N, s0, v0);
        ans = add_phase(ans, v1 + 1, v2, N, s1, v1);
        ans = add_phase(ans, v2 + 1, v3 - 1, N, s2, v2);
    }

    return ans;
}

}  // namespace

int main() {
    assert(G_mod(18) == 530);
    assert(G_mod(1'000'000) == 951'640'919LL);

    std::cout << G_mod(10'000'000'000'000'000ULL) << '\n';
    return 0;
}
