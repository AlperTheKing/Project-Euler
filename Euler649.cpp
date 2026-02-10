#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'000ULL;

u64 mod_mul(u64 a, u64 b) { return static_cast<u64>((u128)a * b % (u128)kMod); }

u64 mod_pow(u64 a, u64 e) {
    u64 r = 1 % kMod;
    while (e) {
        if (e & 1) r = mod_mul(r, a);
        a = mod_mul(a, a);
        e >>= 1;
    }
    return r;
}

std::vector<u64> sg_counts_1d(const int n) {
    std::vector<u64> cnt(8, 0);
    std::vector<std::uint8_t> buf(8, 0);
    const int moves[4] = {2, 3, 5, 7};
    for (int i = 0; i < n; ++i) {
        std::uint8_t g = 0;
        if (i >= 2) {
            bool seen[6] = {false, false, false, false, false, false};
            for (int d : moves) {
                if (i >= d) seen[buf[static_cast<std::size_t>((i - d) & 7)]] = true;
            }
            while (g < 6 && seen[g]) ++g;
        }
        buf[static_cast<std::size_t>(i & 7)] = g;
        ++cnt[g];
    }
    return cnt;
}

u64 M(const int n, const int c) {
    const auto cnt1 = sg_counts_1d(n);
    const int SZ = 8;
    std::vector<u64> cnt2(SZ, 0);
    for (int a = 0; a < SZ; ++a) {
        if (cnt1[a] == 0) continue;
        for (int b = 0; b < SZ; ++b) {
            if (cnt1[b] == 0) continue;
            const int g = a ^ b;
            cnt2[g] += cnt1[a] * cnt1[b];
        }
    }
    std::vector<u64> cnt2m(SZ, 0);
    for (int i = 0; i < SZ; ++i) cnt2m[i] = cnt2[i] % kMod;

    std::vector<u64> dp(SZ, 0), ndp(SZ, 0);
    dp[0] = 1;
    for (int t = 0; t < c; ++t) {
        std::fill(ndp.begin(), ndp.end(), 0);
        for (int x = 0; x < SZ; ++x) {
            const u64 vx = dp[x];
            if (vx == 0) continue;
            for (int g = 0; g < SZ; ++g) {
                ndp[x ^ g] += mod_mul(vx, cnt2m[g]);
                ndp[x ^ g] %= kMod;
            }
        }
        dp.swap(ndp);
    }

    const u64 n2m = mod_mul(static_cast<u64>(n) % kMod, static_cast<u64>(n) % kMod);
    const u64 total = mod_pow(n2m, static_cast<u64>(c));
    return (total + kMod - dp[0]) % kMod;
}

}  // namespace

int main() {
    assert(M(3, 1) == 4);
    assert(M(3, 2) == 40);
    assert(M(9, 3) == 450304);

    const u64 ans = M(10'000'019, 100);
    std::cout << std::setw(9) << std::setfill('0') << ans << "\n";
    return 0;
}

