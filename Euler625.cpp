#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

using u64 = unsigned long long;
using u128 = __uint128_t;

static constexpr int MOD = 998244353;
static constexpr int INV2 = (MOD + 1) / 2;

static inline int mod_add(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static inline int mod_sub(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

static inline int mod_mul(u64 a, u64 b) { return (int)((u128)(a % MOD) * (b % MOD) % MOD); }

static inline int tri_mod(u64 l, u64 r) {
    u64 cnt = r - l + 1;
    int a = (int)((l + r) % MOD);
    int b = (int)(cnt % MOD);
    return mod_mul((u64)a * b, INV2);
}

struct PhiPrefix {
    int limit;
    std::vector<int> pref;
};

static PhiPrefix build_phi_prefix(int limit) {
    std::vector<int> phi(limit + 1, 0);
    std::vector<int> primes;
    std::vector<unsigned char> is_comp(limit + 1, 0);
    primes.reserve((size_t)limit / 10);

    phi[0] = 0;
    if (limit >= 1) phi[1] = 1;

    for (int i = 2; i <= limit; ++i) {
        if (!is_comp[i]) {
            primes.push_back(i);
            phi[i] = i - 1;
        }
        for (int p : primes) {
            long long v = 1LL * p * i;
            if (v > limit) break;
            is_comp[(int)v] = 1;
            if (i % p == 0) {
                phi[(int)v] = phi[i] * p;
                break;
            }
            phi[(int)v] = phi[i] * (p - 1);
        }
    }

    std::vector<int> pref(limit + 1, 0);
    for (int i = 1; i <= limit; ++i) {
        pref[i] = mod_add(pref[i - 1], phi[i] % MOD);
    }
    return PhiPrefix{limit, std::move(pref)};
}

static int sum_phi(u64 n, const PhiPrefix &pre, std::unordered_map<u64, int> &memo) {
    if (n <= (u64)pre.limit) return pre.pref[(size_t)n];
    auto it = memo.find(n);
    if (it != memo.end()) return it->second;

    int res = mod_mul(n, n + 1);
    res = mod_mul(res, INV2);

    for (u64 l = 2; l <= n;) {
        u64 q = n / l;
        u64 r = n / q;
        int cnt = (int)((r - l + 1) % MOD);
        int sub = mod_mul(cnt, (u64)sum_phi(q, pre, memo));
        res = mod_sub(res, sub);
        l = r + 1;
    }

    memo.emplace(n, res);
    return res;
}

static int G(u64 N, const PhiPrefix &pre, std::unordered_map<u64, int> &memo) {
    (void)sum_phi(N, pre, memo);
    int ans = 0;
    for (u64 l = 1; l <= N;) {
        u64 q = N / l;
        u64 r = N / q;
        int sum_d = tri_mod(l, r);
        int sphi = sum_phi(q, pre, memo);
        ans = mod_add(ans, mod_mul((u64)sum_d, (u64)sphi));
        l = r + 1;
    }
    return ans;
}

static u64 G_bruteforce(u64 N) {
    u64 ans = 0;
    for (u64 j = 1; j <= N; ++j) {
        for (u64 i = 1; i <= j; ++i) ans += std::gcd(i, j);
    }
    return ans;
}

int main() {
    const int LIMIT = 5'000'000;
    const PhiPrefix pre = build_phi_prefix(LIMIT);
    std::unordered_map<u64, int> memo;
    memo.reserve(1 << 18);

    assert(G_bruteforce(10) == 122);
    assert(G(10, pre, memo) == 122);
    std::cout << G(100'000'000'000ULL, pre, memo) << "\n";
    return 0;
}
