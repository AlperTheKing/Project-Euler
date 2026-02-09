#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>

using namespace std;

using u64 = uint64_t;
using u128 = unsigned __int128;

static constexpr u64 MOD = 987654321ULL;

static inline u64 mod_add(u64 a, u64 b) {
    u64 s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static inline u64 mod_sub(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + MOD - b);
}

static inline u64 mod_mul(u64 a, u64 b) {
    return (u64)((u128)(a % MOD) * (u128)(b % MOD) % MOD);
}

struct Solver {
    unordered_map<u64, u64> memoP;
    unordered_map<u64, u64> memoS;

    u64 P(u64 n) {
        if (n <= 1) return 1;
        auto it = memoP.find(n);
        if (it != memoP.end()) return it->second;
        u64 m = n / 2;
        u64 pm = P(m);
        // P(n) = 2 * (1 + floor(n/2) - P(floor(n/2)))
        u64 res = 2 * (1 + m - pm);
        memoP.emplace(n, res);
        return res;
    }

    u64 S(u64 n) {
        if (n == 0) return 0;
        if (n == 1) return 1 % MOD;
        auto it = memoS.find(n);
        if (it != memoS.end()) return it->second;
        u64 m = n / 2;
        u64 Sm = S(m);

        u64 res;
        if (n & 1ULL) {
            // n = 2m+1:
            // S(2m+1) = 1 + 2m(m+3) - 4 S(m)
            u64 term = mod_mul(2 % MOD, mod_mul(m % MOD, (m + 3) % MOD));
            u64 sub4 = mod_mul(4 % MOD, Sm);
            res = mod_add(1 % MOD, mod_sub(term, sub4));
        } else {
            // n = 2m:
            // S(2m) = 2m^2 + 4m - 1 - 4 S(m) + 2 P(m)
            u64 mmod = m % MOD;
            u64 term2m2 = mod_mul(2 % MOD, mod_mul(mmod, mmod));
            u64 term4m = mod_mul(4 % MOD, mmod);
            u64 base = mod_sub(mod_add(term2m2, term4m), 1 % MOD);
            u64 sub4 = mod_mul(4 % MOD, Sm);
            u64 add2P = mod_mul(2 % MOD, P(m) % MOD);
            res = mod_add(mod_sub(base, sub4), add2P);
        }

        memoS.emplace(n, res);
        return res;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Solver solver;

    // Validation points from the problem statement.
    assert(solver.P(1) == 1);
    assert(solver.P(9) == 6);
    assert(solver.P(1000) == 510);
    assert(solver.S(1000) == 268271ULL);

    u64 ans = solver.S(1000000000000000000ULL) % MOD;
    cout << ans << '\n';
    return 0;
}
