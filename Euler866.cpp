#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;
using u128 = unsigned __int128;

static constexpr i64 MOD = 987'654'319LL;

static i64 mod_mul(i64 a, i64 b) {
    return static_cast<i64>((static_cast<u128>(a) * b) % MOD);
}

static i64 mod_inv(i64 a) {
    i64 r0 = a, r1 = MOD;
    i64 s0 = 1, s1 = 0;
    while (r1 != 0) {
        i64 q = r0 / r1;
        i64 r2 = r0 - q * r1;
        r0 = r1;
        r1 = r2;
        i64 s2 = s0 - q * s1;
        s0 = s1;
        s1 = s2;
    }
    if (s0 < 0) s0 += MOD;
    return s0;
}

static i64 G_mod(int N) {
    std::vector<i64> g(N + 1, 0);
    g[0] = 1;

    for (int n = 1; n <= N; ++n) {
        i64 sum = 0;
        for (int i = 1; i <= n; ++i) {
            sum += mod_mul(g[i - 1], g[n - i]);
            if (sum >= MOD) sum -= MOD;
        }

        i64 H = static_cast<i64>(n) * (2LL * n - 1) % MOD;
        i64 invn = mod_inv(n % MOD);
        g[n] = mod_mul(mod_mul(H, sum), invn);
    }

    return g[N];
}

int main() {
    assert(G_mod(4) == 994);
    std::cout << G_mod(100) << '\n';
    return 0;
}
