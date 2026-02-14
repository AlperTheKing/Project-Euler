#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <functional>

using namespace std;

using u64 = uint64_t;
using u128 = unsigned __int128;

static u64 isqrt_u64(u64 x) {
    long double r = sqrt((long double)x);
    u64 y = (u64)r;
    while ((u128)(y + 1) * (y + 1) <= (u128)x) ++y;
    while ((u128)y * y > (u128)x) --y;
    return y;
}

static u64 brute_P(u64 N) {
    // Brute for validation (N up to 1e6 is fine).
    u64 cnt = 0;
    u64 mmax = isqrt_u64(N - 1);
    for (u64 m = 2; m <= mmax; m++) {
        u64 mm = m * m;
        u64 nmax = isqrt_u64(N - mm);
        if (nmax >= m) nmax = m - 1;
        for (u64 n = 1; n <= nmax; n++) {
            if (((m ^ n) & 1ULL) == 0) continue; // opposite parity
            if (std::gcd(m, n) != 1) continue;
            cnt++;
        }
    }
    return cnt;
}

struct LinearSieve {
    int n;
    vector<uint32_t> spf;
    vector<int> primes;

    explicit LinearSieve(int n_) : n(n_), spf(n_ + 1, 0) {
        primes.reserve((size_t)(n / 10));
        for (int i = 2; i <= n; i++) {
            if (spf[i] == 0) {
                spf[i] = (uint32_t)i;
                primes.push_back(i);
            }
            for (int p : primes) {
                long long v = 1LL * p * i;
                if (v > n) break;
                spf[(int)v] = (uint32_t)p;
                if (p == (int)spf[i]) break;
            }
        }
    }

    inline void distinct_prime_factors(int x, int *buf, int &k) const {
        k = 0;
        while (x > 1) {
            int p = (int)spf[x];
            buf[k++] = p;
            while (x % p == 0) x /= p;
        }
    }
};

static inline u64 coprime_count(u64 t, int x, const LinearSieve &sv, int *pf_buf) {
    // Count 1<=n<=t with gcd(n,x)=1 via inclusion-exclusion.
    if (t == 0) return 0;
    if (x == 1) return t;

    int k = 0;
    sv.distinct_prime_factors(x, pf_buf, k);

    u64 res = 0;
    int masks = 1 << k;
    for (int mask = 0; mask < masks; mask++) {
        u64 prod = 1;
        int bits = 0;
        for (int i = 0; i < k; i++) {
            if (mask & (1 << i)) {
                prod *= (u64)pf_buf[i];
                bits++;
            }
        }
        u64 term = t / prod;
        if ((bits & 1) == 0) res += term;
        else res -= term;
    }
    return res;
}

static u64 P_large(u64 N) {
    // Count primitive Pythagorean triples with a<b<c<=N.
    // Euclid: c = m^2 + n^2, m>n, gcd(m,n)=1, m-n odd.

    u64 mmax = isqrt_u64(N - 1);
    LinearSieve sv((int)mmax);

    int pf_buf[16];

    u64 total = 0;
    u64 nmax = mmax;

    for (u64 m = 2; m <= mmax; m++) {
        u128 mm = (u128)m * (u128)m;
        while (nmax >= 1 && mm + (u128)nmax * (u128)nmax > (u128)N) --nmax;
        if (nmax == 0) break;

        u64 t = nmax;
        if (t > m - 1) t = m - 1;
        if (t == 0) continue;

        if (m & 1ULL) {
            // m odd -> n even: n=2x, x<=t/2, gcd(m,x)=1
            total += coprime_count(t / 2, (int)m, sv, pf_buf);
        } else {
            // m even -> n odd. Let q=m/2, need gcd(q,n)=1 and n odd.
            int q = (int)(m >> 1);
            u64 all = coprime_count(t, q, sv, pf_buf);
            if ((q & 1) == 0) {
                // q even => all coprime n are odd.
                total += all;
            } else {
                // q odd: subtract even n=2x with gcd(q,x)=1.
                total += all - coprime_count(t / 2, q, sv, pf_buf);
            }
        }
    }

    return total;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Validation points from the problem statement.
    assert(brute_P(20) == 3);
    assert(brute_P(1'000'000ULL) == 159139ULL);

    const u64 N = 3141592653589793ULL;
    cout << P_large(N) << '\n';
    return 0;
}
