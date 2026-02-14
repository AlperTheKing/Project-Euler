#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>

// Project Euler 606: only numbers of the form (p*q)^3 (p<q primes) have 252 gozinta chains,
// so S(N) reduces to sum_{p<q, pq <= cbrt(N)} (pq)^3. For cbrt(10^36)=10^12 we compute this
// modulo 10^9 using a Min_25 style prime power-sum sieve for sum_{p<=x} p^3.

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

static constexpr u64 MOD = 1000000000ULL;  // last 9 digits

static inline u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}
static inline u64 sub_mod(u64 a, u64 b) { return (a >= b) ? (a - b) : (a + MOD - b); }
static inline u64 mul_mod(u64 a, u64 b) { return (u64)((u128)a * b % MOD); }

static u64 sum_cubes_mod(u64 n) {
    // sum_{i=1..n} i^3 = (n(n+1)/2)^2
    const u128 t = (u128)n * (n + 1) / 2;
    const u64 tm = (u64)(t % MOD);
    return mul_mod(tm, tm);
}

static std::vector<u64> sieve_primes(u64 n) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (u64 p = 2; p * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (u64 j = p * p; j <= n; j += p) is_prime[j] = false;
    }
    std::vector<u64> primes;
    primes.reserve(n / 10);
    for (u64 i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

struct PrimeCubeSummatory {
    u64 n = 0;
    u64 sq = 0;
    std::vector<u64> vals;   // distinct floor(n/i)
    std::vector<u64> g;      // sum_{p<=vals[j]} p^3 (mod MOD) after sieve
    std::vector<int> id1;    // vals index for x<=sq
    std::vector<int> id2;    // vals index for x>sq, keyed by n/x
    std::vector<u64> primes; // primes up to sq
    std::vector<u64> pref_p3;

    explicit PrimeCubeSummatory(u64 n_) : n(n_) {
        sq = (u64)std::sqrt((long double)n);
        primes = sieve_primes(sq);

        for (u64 l = 1; l <= n;) {
            const u64 w = n / l;
            vals.push_back(w);
            l = n / w + 1;
        }

        const int m = (int)vals.size();
        g.resize(m);
        id1.assign((size_t)sq + 1, -1);
        id2.assign((size_t)sq + 1, -1);
        for (int i = 0; i < m; ++i) {
            const u64 w = vals[i];
            g[i] = sub_mod(sum_cubes_mod(w), 1); // exclude 1^3
            if (w <= sq) id1[w] = i;
            else id2[n / w] = i;
        }

        pref_p3.assign(primes.size() + 1, 0);
        for (std::size_t i = 0; i < primes.size(); ++i) {
            const u64 p = primes[i];
            const u64 p3 = (u64)((u128)p * p % MOD * p % MOD);
            pref_p3[i + 1] = add_mod(pref_p3[i], p3);
        }

        for (std::size_t i = 0; i < primes.size(); ++i) {
            const u64 p = primes[i];
            const u64 p2 = p * p;
            if (p2 > n) break;
            const u64 p3 = (u64)((u128)p * p % MOD * p % MOD);

            for (int j = 0; j < m && vals[j] >= p2; ++j) {
                const u64 w = vals[j];
                const int idx = id(w / p);
                const u64 term = sub_mod(g[idx], pref_p3[i]);
                g[j] = sub_mod(g[j], mul_mod(p3, term));
            }
        }
    }

    inline int id(u64 x) const { return (x <= sq) ? id1[x] : id2[n / x]; }
    inline u64 sum_p3(u64 x) const { return g[id(x)]; }
};

static u64 g_from_pattern(const std::vector<int>& exps, u64 limit) {
    if (exps.empty()) return 1;
    int omega = 0;
    for (int e : exps) omega += e;

    std::vector<std::vector<u64>> C((size_t)omega + 1, std::vector<u64>((size_t)omega + 1, 0));
    C[0][0] = 1;
    for (int n = 1; n <= omega; ++n) {
        C[n][0] = C[n][n] = 1;
        for (int k = 1; k < n; ++k) C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
    }

    auto binom_u64 = [](int n, int k) -> u64 {
        if (k < 0 || k > n) return 0;
        k = std::min(k, n - k);
        u128 r = 1;
        for (int i = 1; i <= k; ++i) r = (r * (u128)(n - k + i)) / (u128)i;
        return (u64)r;
    };

    std::vector<u128> total((size_t)omega + 1, 0);
    for (int k = 1; k <= omega; ++k) {
        u128 prod = 1;
        for (int a : exps) prod *= (u128)binom_u64(a + k - 1, k - 1);
        total[(size_t)k] = prod;
    }

    u128 g = 0;
    for (int k = 1; k <= omega; ++k) {
        i128 fk = 0;
        for (int i = 1; i <= k; ++i) {
            const u128 term = (u128)C[k][i] * total[(size_t)i];
            if (((k - i) & 1) != 0) fk -= (i128)term;
            else fk += (i128)term;
        }
        g += (u128)fk;
        if (g > (u128)limit) return limit + 1;
    }
    return (u64)g;
}

static void gen_partitions(int rem, int max_part, std::vector<int>& cur, std::vector<std::vector<int>>& out) {
    if (rem == 0) {
        out.push_back(cur);
        return;
    }
    for (int x = std::min(rem, max_part); x >= 1; --x) {
        cur.push_back(x);
        gen_partitions(rem - x, x, cur, out);
        cur.pop_back();
    }
}

static u128 brute_sum_pairs_cubed(u64 X) {
    std::vector<u64> primes = sieve_primes(X);
    u128 sum = 0;
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const u64 p = primes[i];
        for (std::size_t j = i + 1; j < primes.size(); ++j) {
            const u64 q = primes[j];
            if ((u128)p * q > X) break;
            const u128 pq = (u128)p * q;
            sum += pq * pq * pq;
        }
    }
    return sum;
}

int main() {
    // Verify the exponent-pattern characterization of 252.
    std::vector<std::vector<int>> patterns;
    for (int omega = 1; omega <= 20; ++omega) {
        std::vector<std::vector<int>> parts;
        std::vector<int> cur;
        gen_partitions(omega, omega, cur, parts);
        for (auto& p : parts) {
            std::sort(p.begin(), p.end(), std::greater<int>());
            if (g_from_pattern(p, 252) == 252) patterns.push_back(p);
        }
    }
    std::sort(patterns.begin(), patterns.end());
    patterns.erase(std::unique(patterns.begin(), patterns.end()), patterns.end());
    assert(patterns.size() == 1 && patterns[0].size() == 2 && patterns[0][0] == 3 && patterns[0][1] == 3);

    // Statement validations.
    assert((u64)brute_sum_pairs_cubed(100) == 8462952ULL); // S(10^6)
    assert((u64)brute_sum_pairs_cubed(10000) == 623291998881978ULL); // S(10^12)

    const u64 X = 1000000000000ULL; // cbrt(10^36)
    PrimeCubeSummatory pc(X);

    u64 ans = 0;
    for (u64 p : pc.primes) { // primes up to sqrt(X)
        const u64 lim = X / p;
        if (lim <= p) break;
        const u64 sum_q = sub_mod(pc.sum_p3(lim), pc.sum_p3(p));
        const u64 p3 = (u64)((u128)p * p % MOD * p % MOD);
        ans = add_mod(ans, mul_mod(p3, sum_q));
    }

    std::cout << std::setw(9) << std::setfill('0') << ans << "\n";
    return 0;
}
