#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_add(u64 a, u64 b) {
    a += b;
    if (a >= kMod) a -= kMod;
    return a;
}

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

struct Sieve {
    int n;
    std::vector<int> primes;
    std::vector<int> spf;
};

Sieve sieve_spf(const int n) {
    Sieve s;
    s.n = n;
    s.spf.assign(static_cast<std::size_t>(n + 1), 0);
    for (int i = 2; i <= n; ++i) {
        if (s.spf[i] == 0) {
            s.spf[i] = i;
            s.primes.push_back(i);
        }
        for (int p : s.primes) {
            if (p > s.spf[i] || (u64)i * (u64)p > (u64)n) break;
            s.spf[i * p] = p;
        }
    }
    return s;
}

u64 solve(const int N) {
    const Sieve sv = sieve_spf(N);
    const int P = static_cast<int>(sv.primes.size());

    std::vector<int> pidx(static_cast<std::size_t>(N + 1), -1);
    for (int i = 0; i < P; ++i) pidx[static_cast<std::size_t>(sv.primes[i])] = i;

    std::vector<u64> invp(P, 0), invpm1(P, 0);
    for (int i = 0; i < P; ++i) {
        const u64 p = static_cast<u64>(sv.primes[i]);
        invp[i] = mod_pow(p % kMod, kMod - 2);
        invpm1[i] = mod_pow((p - 1) % kMod, kMod - 2);
    }

    std::vector<u64> invPowA(P, 1);  // p^{-v_p((n)!)} after processing n
    std::vector<u64> powE1(P, 0);    // p^{e_p(B(n))+1}
    for (int i = 0; i < P; ++i) powE1[i] = static_cast<u64>(sv.primes[i]) % kMod;

    u64 S = 0;

    for (int n = 1; n <= N; ++n) {
        for (int i = 0; i < P; ++i) powE1[i] = mod_mul(powE1[i], invPowA[i]);

        int x = n;
        while (x > 1) {
            const int p = sv.spf[x];
            int v = 0;
            while (x % p == 0) {
                x /= p;
                ++v;
            }
            const int i = pidx[static_cast<std::size_t>(p)];
            const u64 add_e = static_cast<u64>(n - 1) * static_cast<u64>(v);
            if (add_e) powE1[i] = mod_mul(powE1[i], mod_pow(static_cast<u64>(p) % kMod, add_e));
            invPowA[i] = mod_mul(invPowA[i], mod_pow(invp[i], static_cast<u64>(v)));
        }

        u64 D = 1;
        for (int i = 0; i < P; ++i) {
            const u64 t = powE1[i];
            const u64 num = (t + kMod - 1) % kMod;
            const u64 term = mod_mul(num, invpm1[i]);
            D = mod_mul(D, term);
        }

        S = mod_add(S, D);
    }

    return S;
}

}  // namespace

int main() {
    assert(solve(5) == 5736);
    assert(solve(10) == (141740594713218418ULL % kMod));
    assert(solve(100) == 332792866);
    std::cout << solve(20'000) << "\n";
    return 0;
}
