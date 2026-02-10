#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

u64 floor_sqrt_u128(const u128 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((u128)(r + 1) * (r + 1) <= n) ++r;
    while ((u128)r * r > n) --r;
    return r;
}

i128 Fk(const int k, const u64 N) {
    if ((k & 1) == 0 || k < 3) return 0;
    const u64 d = static_cast<u64>(k - 2);
    const u64 c = static_cast<u64>(k - 4) * static_cast<u64>(k - 4);

    const u64 sqrtN = isqrt_u64(N);
    if (sqrtN <= 1) return 0;

    const u64 tA = (sqrtN - 1) / (2 * d);

    const u64 M = (2 * N) / c;
    const u128 disc = 1 + (u128)4 * d * (u128)M;
    const u64 s = floor_sqrt_u128(disc);
    const u64 tB = (s > 0) ? (static_cast<u64>((s - 1) / (2 * d))) : 0;

    const u64 T = std::min(tA, tB);
    if (T == 0) return 0;

    const i128 TT = static_cast<i128>(T);
    const i128 S1 = TT * (TT + 1) / 2;
    const i128 S2 = TT * (TT + 1) * (2 * TT + 1) / 6;

    const i128 dd = static_cast<i128>(d);
    const i128 cc = static_cast<i128>(c);

    const i128 sumA = TT + 4 * dd * S1 + 4 * dd * dd * S2;
    const i128 sum_dt2_t = dd * S2 + S1;
    const i128 sumB = cc * (sum_dt2_t / 2);
    return sumA + sumB;
}

i128 sum_all(const u64 N) {
    i128 ans = 0;
    for (int k = 3;; k += 2) {
        const u64 d = static_cast<u64>(k - 2);
        const u64 c = static_cast<u64>(k - 4) * static_cast<u64>(k - 4);
        const u128 A1 = (u128)(2ULL * (u64)k - 3) * (2ULL * (u64)k - 3);
        const u128 B1 = (u128)c * (u64)(k - 1) / 2;
        if (A1 > N || B1 > N) break;
        ans += Fk(k, N);
    }
    return ans;
}

std::string to_string_i128(i128 x) {
    if (x == 0) return "0";
    bool neg = false;
    if (x < 0) {
        neg = true;
        x = -x;
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

}  // namespace

int main() {
    assert(Fk(3, 100) == 184);
    assert(sum_all(1000) == 14993);

    std::cout << to_string_i128(sum_all(1'000'000'000'000ULL)) << "\n";
    return 0;
}

