#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 MOD = 998244353ULL;

static u64 mod_pow(u64 a, u64 e) {
    u64 r = 1;
    while (e > 0) {
        if (e & 1ULL) r = (static_cast<u128>(r) * a) % MOD;
        a = (static_cast<u128>(a) * a) % MOD;
        e >>= 1ULL;
    }
    return r;
}

static u64 T_mod(int N) {
    std::vector<u64> U(N + 1, 0), V(N + 1, 0), E(N + 1, 0);
    E[0] = 1;

    for (int n = 1; n <= N; ++n) {
        u64 sum_uv = 0;
        u64 sum_uu = 0;
        for (int i = 0; i < n; ++i) {
            sum_uv = (sum_uv + static_cast<u128>(U[i]) * V[n - 1 - i]) % MOD;
            sum_uu = (sum_uu + static_cast<u128>(U[i]) * U[n - 1 - i]) % MOD;
        }

        u64 vn = (n == 1 ? 1ULL : 0ULL);
        vn = (vn + 9ULL * sum_uv) % MOD;
        V[n] = vn;

        u64 un = (V[n - 1] + 9ULL * sum_uu) % MOD;
        U[n] = un;
    }

    for (int n = 1; n <= N; ++n) {
        u64 sum_ue = 0;
        for (int i = 0; i < n; ++i) {
            sum_ue = (sum_ue + static_cast<u128>(U[i]) * E[n - 1 - i]) % MOD;
        }
        E[n] = (10ULL * sum_ue) % MOD;
    }

    u64 sumE = 0;
    for (int n = 1; n <= N; ++n) {
        sumE += E[n];
        if (sumE >= MOD) sumE -= MOD;
    }

    u64 inv10 = mod_pow(10, MOD - 2);
    return static_cast<u64>((static_cast<u128>(9) * inv10 % MOD) * sumE % MOD);
}

static u64 T_exact_small(int N) {
    std::vector<u128> U(N + 1, 0), V(N + 1, 0), E(N + 1, 0);
    E[0] = 1;

    for (int n = 1; n <= N; ++n) {
        u128 sum_uv = 0;
        u128 sum_uu = 0;
        for (int i = 0; i < n; ++i) {
            sum_uv += U[i] * V[n - 1 - i];
            sum_uu += U[i] * U[n - 1 - i];
        }

        V[n] = (n == 1 ? 1 : 0) + 9 * sum_uv;
        U[n] = V[n - 1] + 9 * sum_uu;
    }

    for (int n = 1; n <= N; ++n) {
        u128 sum_ue = 0;
        for (int i = 0; i < n; ++i) sum_ue += U[i] * E[n - 1 - i];
        E[n] = 10 * sum_ue;
    }

    u128 ans = 0;
    for (int n = 1; n <= N; ++n) {
        ans += (9 * E[n]) / 10;
    }
    return static_cast<u64>(ans);
}

int main() {
    assert(T_exact_small(6) == 261ULL);
    assert(T_exact_small(30) == 5'576'195'181'577'716ULL);

    std::cout << T_mod(10'000) << '\n';
    return 0;
}
