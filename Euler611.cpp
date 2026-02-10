#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 611: count n <= N where #{0<a<b : a^2 + b^2 = n} is odd.
// Reduce to two arithmetic families via r2(n) and prime-exponent parity; count with a prime-counting sieve for p ≡ 1 (mod 4).

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

static u64 isqrt_u64(u64 x) {
    u64 r = (u64)std::sqrt((long double)x);
    while ((u128)(r + 1) * (r + 1) <= x) ++r;
    while ((u128)r * r > x) --r;
    return r;
}

static int ilog2_u64(u64 x) {
    int r = 0;
    while ((1ULL << (r + 1)) <= x) ++r;
    return r;
}

static std::vector<int> sieve_primes_int(int n) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; (u64)p * p <= (u64)n; ++p) {
        if (!is_prime[p]) continue;
        for (u64 j = (u64)p * p; j <= (u64)n; j += (u64)p) is_prime[(size_t)j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

struct PrimeCountMod4 {
    u64 n = 0;
    u64 sq = 0;
    std::vector<u64> vals;
    std::vector<u64> g_pi;
    std::vector<i64> g_chi;
    std::vector<int> id1;
    std::vector<int> id2;
    std::vector<int> primes;
    std::vector<i64> pref_chi_prime;

    static i64 chi_int(u64 x) {
        if ((x & 1ULL) == 0) return 0;
        return (x & 3ULL) == 1 ? 1 : -1;
    }

    static i64 sum_chi_1_to(u64 m) {
        // chi(k)=0 for even; for odd k, +1 if 1 mod4 else -1.
        const u64 c1 = (m + 3) / 4;
        const u64 c3 = (m + 1) / 4;
        return (i64)c1 - (i64)c3;
    }

    explicit PrimeCountMod4(u64 n_) : n(n_) {
        sq = (u64)std::sqrt((long double)n);
        primes = sieve_primes_int((int)sq);

        for (u64 l = 1; l <= n;) {
            const u64 w = n / l;
            vals.push_back(w);
            l = n / w + 1;
        }
        const int m = (int)vals.size();

        g_pi.resize(m);
        g_chi.resize(m);
        id1.assign((size_t)sq + 1, -1);
        id2.assign((size_t)sq + 1, -1);

        for (int i = 0; i < m; ++i) {
            const u64 w = vals[i];
            g_pi[i] = (w >= 2) ? (w - 1) : 0; // count of integers in [2..w]
            i64 s = (w >= 1) ? (sum_chi_1_to(w) - 1) : 0; // sum_{k=2..w} chi(k)
            g_chi[i] = s;
            if (w <= sq) id1[w] = i;
            else id2[n / w] = i;
        }

        pref_chi_prime.assign(primes.size() + 1, 0);
        for (std::size_t i = 0; i < primes.size(); ++i) {
            const int p = primes[i];
            pref_chi_prime[i + 1] = pref_chi_prime[i] + chi_int((u64)p);
        }

        for (std::size_t i = 0; i < primes.size(); ++i) {
            const u64 p = (u64)primes[i];
            const u64 p2 = p * p;
            if (p2 > n) break;
            const i64 chi_p = chi_int(p);

            for (int j = 0; j < m && vals[j] >= p2; ++j) {
                const u64 w = vals[j];
                const int idx = id(w / p);
                g_pi[j] -= g_pi[idx] - (u64)i;
                if (chi_p != 0) g_chi[j] -= chi_p * (g_chi[idx] - pref_chi_prime[i]);
            }
        }
    }

    inline int id(u64 x) const { return (x <= sq) ? id1[x] : id2[n / x]; }

    inline u64 pi(u64 x) const { return g_pi[id(x)]; }

    inline i64 sum_chi_primes(u64 x) const { return g_chi[id(x)]; }

    inline u64 pi1(u64 x) const {
        if (x < 5) return 0;
        const u64 pix = pi(x);
        const u64 odd_primes = pix - 1; // exclude prime 2
        const i64 s = sum_chi_primes(x);
        return (u64)((odd_primes + s) / 2);
    }
};

static std::vector<int> sieve_spf(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const long long v = 1LL * p * i;
            if (v > n) break;
            spf[(int)v] = p;
            if (p == spf[i]) break;
        }
    }
    return spf;
}

static u64 count_B(u64 N, const std::vector<uint8_t>& parity_odd1mod4) {
    const u64 lim = isqrt_u64(N);
    u64 ans = 0;
    // Unique parameterization for squares / twice-squares: n = 2^a * m^2 with m odd.
    for (u64 m = 1; m <= lim; m += 2) {
        if (!parity_odd1mod4[(size_t)m]) continue;
        const u64 base = m * m;
        const u64 t = N / base;
        ans += (u64)(ilog2_u64(t) + 1);
    }
    return ans;
}

static u64 count_A_even_v2(u64 N, const PrimeCountMod4& pc, const std::vector<int>& primes_small) {
    // a=0 term: sum_{p≡1(4)} (T - T/p), T=floor(sqrt(N/p)).
    u128 S1 = 0;
    const u64 tmax = isqrt_u64(N);
    for (u64 t = 1; t <= tmax; ++t) {
        const u64 R = N / (t * t);
        if (R < 5) break;
        const u64 L = N / ((t + 1) * (t + 1));
        const u64 cnt = pc.pi1(R) - pc.pi1(L);
        S1 += (u128)t * cnt;
    }

    u128 corr = 0;
    const u64 cbrtN = 10000; // for N up to 1e12 (exact)
    for (int p : primes_small) {
        if ((u64)p > cbrtN) break;
        if (p % 4 != 1) continue;
        const u64 T = isqrt_u64(N / (u64)p);
        corr += (u64)(T / (u64)p);
    }

    u128 total = S1 - corr;

    // higher exponents: p^{4a+1} with a>=1 (start at p^5)
    for (int p : primes_small) {
        if (p % 4 != 1) continue;
        u128 pow = 1;
        for (int i = 0; i < 5; ++i) pow *= (u64)p;
        if (pow > N) break;
        const u128 p4 = (u128)p * p * p * p;
        while (pow <= N) {
            const u64 T = isqrt_u64((u64)(N / (u64)pow));
            total += (u128)(T - T / (u64)p);
            if (pow > (u128)N / p4) break;
            pow *= p4;
        }
    }

    return (u64)total;
}

static u64 F(u64 N, const PrimeCountMod4& pc, const std::vector<int>& primes_small, const std::vector<int>& spf,
             const std::vector<uint8_t>& parity_odd1mod4) {
    (void)spf;
    const u64 a_even = count_A_even_v2(N, pc, primes_small);
    const u64 a_odd = count_A_even_v2(N / 2, pc, primes_small); // multiply-by-2 bijection
    return a_even + a_odd + count_B(N, parity_odd1mod4);
}

int main() {
    const u64 Nmax = 1000000000000ULL;
    PrimeCountMod4 pc(Nmax);

    const int LIM = 1000000;
    const std::vector<int> spf = sieve_spf(LIM);
    const std::vector<int> primes_small = sieve_primes_int(LIM);

    // parity[m] = 1 iff m has an odd number of primes ≡1 (mod 4) to odd exponent.
    std::vector<uint8_t> parity((size_t)LIM + 1, 0);
    parity[1] = 0;
    for (int x = 2; x <= LIM; ++x) {
        int n = x;
        uint8_t par = 0;
        while (n > 1) {
            const int p = spf[n];
            int e = 0;
            while (n % p == 0) {
                n /= p;
                e ^= 1;
            }
            if (e && (p % 4 == 1)) par ^= 1;
        }
        parity[x] = par;
    }

    // Statement validations.
    assert(F(5, pc, primes_small, spf, parity) == 1ULL);
    assert(F(100, pc, primes_small, spf, parity) == 27ULL);
    assert(F(1000, pc, primes_small, spf, parity) == 233ULL);
    assert(F(1000000ULL, pc, primes_small, spf, parity) == 112168ULL);

    std::cout << F(Nmax, pc, primes_small, spf, parity) << "\n";
    return 0;
}
