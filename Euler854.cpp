#include <cassert>
#include <cstdint>
#include <iostream>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 MOD = 1'234'567'891ULL;

static u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % MOD);
}

static u64 M_small(int p) {
    if (p == 3) return 2;
    if (p == 2) return 1;
    if (p & 1) return 1;

    int k = p / 2;
    u64 f0 = 0, f1 = 1;
    u64 l0 = 2, l1 = 1;
    for (int i = 2; i <= k; ++i) {
        u64 fn = f0 + f1;
        f0 = f1;
        f1 = fn;
        u64 ln = l0 + l1;
        l0 = l1;
        l1 = ln;
    }
    u64 fk = (k == 1 ? 1 : f1);
    u64 lk = (k == 1 ? 1 : l1);
    return (k & 1) ? lk : fk;
}

static u64 compute_P(int n) {
    u64 ans = 1;
    if (n >= 3) ans = mul_mod(ans, 2);

    int half = n / 2;
    u64 f_prev = 0, f_cur = 1;
    u64 l_prev = 2, l_cur = 1;

    for (int k = 1; k <= half; ++k) {
        u64 fk = (k == 1 ? 1 : 0);
        u64 lk = (k == 1 ? 1 : 0);

        if (k >= 2) {
            u64 fn = (f_prev + f_cur) % MOD;
            f_prev = f_cur;
            f_cur = fn;
            fk = f_cur;

            u64 ln = (l_prev + l_cur) % MOD;
            l_prev = l_cur;
            l_cur = ln;
            lk = l_cur;
        }

        if (k >= 2) {
            ans = mul_mod(ans, (k & 1) ? lk : fk);
        }
    }

    return ans;
}

int main() {
    assert(M_small(18) == 76);
    assert(compute_P(10) == 264);

    std::cout << compute_P(1'000'000) << '\n';
    return 0;
}
