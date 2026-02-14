#include <cstdint>
#include <iostream>
#include <unordered_map>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    if (mod == 1) return 0;
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) r = static_cast<u64>((static_cast<u128>(r) * a) % mod);
        a = static_cast<u64>((static_cast<u128>(a) * a) % mod);
        e >>= 1ULL;
    }
    return r;
}

static u64 inv_mod(u64 a, u64 mod) {
    i64 t = 0, nt = 1;
    i64 r = static_cast<i64>(mod), nr = static_cast<i64>(a % mod);
    while (nr != 0) {
        const i64 q = r / nr;
        const i64 tt = t - q * nt;
        t = nt;
        nt = tt;
        const i64 rr = r - q * nr;
        r = nr;
        nr = rr;
    }
    if (t < 0) t += static_cast<i64>(mod);
    return static_cast<u64>(t);
}

static u64 pow_mod_signed_2(i64 exp, u64 mod) {
    if (mod == 1) return 0;
    if (exp >= 0) return pow_mod(2, static_cast<u64>(exp), mod);
    const u64 inv2 = inv_mod(2, mod);
    return pow_mod(inv2, static_cast<u64>(-exp), mod);
}

static u64 phi(u64 n) {
    if (n == 0) return 0;
    u64 result = n;
    if ((n & 1ULL) == 0) {
        result -= result / 2;
        while ((n & 1ULL) == 0) n >>= 1ULL;
    }
    for (u64 p = 3; p * p <= n; p += 2) {
        if (n % p != 0) continue;
        result -= result / p;
        while (n % p == 0) n /= p;
    }
    if (n > 1) result -= result / n;
    return result;
}

static u64 tower2_fixpoint_mod(u64 mod, std::unordered_map<u64, u64>& memo) {
    auto it = memo.find(mod);
    if (it != memo.end()) return it->second;

    u64 odd = mod;
    i64 twos = 0;
    while ((odd & 1ULL) == 0) {
        odd >>= 1ULL;
        ++twos;
    }

    const u64 t = phi(odd);
    const i64 w = (t == 1) ? 0 : static_cast<i64>(tower2_fixpoint_mod(t, memo)) - twos;
    const u64 odd_part = pow_mod_signed_2(w, odd);
    const u64 ans = (odd_part << twos) % mod;
    memo.emplace(mod, ans);
    return ans;
}

int main() {
    const u64 MOD = 1'000'000'000'000'000ULL;
    std::unordered_map<u64, u64> memo;
    memo.reserve(64);
    memo.emplace(1, 0);
    const u64 top = tower2_fixpoint_mod(MOD, memo);
    std::cout << (top + MOD - 3) % MOD << '\n';
    return 0;
}
