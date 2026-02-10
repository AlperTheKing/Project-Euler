#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 609: group u0 by m=pi(u0) (an interval with exactly one prime), and for each m
// walk the iterated pi-chain m,pi(m),...,1 to count prefix nonprimes.

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = unsigned __int128;

static constexpr u64 MOD = 1000000007ULL;

static std::vector<bool> sieve_is_prime(int n, std::vector<int>& primes) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; (u64)p * p <= (u64)n; ++p) {
        if (!is_prime[p]) continue;
        for (u64 j = (u64)p * p; j <= (u64)n; j += (u64)p) is_prime[(size_t)j] = false;
    }
    primes.clear();
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return is_prime;
}

static u64 P(u64 n) {
    std::vector<int> primes;
    const std::vector<bool> is_prime = sieve_is_prime((int)n, primes);
    const u64 M = (u64)primes.size(); // pi(n)

    // pi(x) needed only for x<=M.
    std::vector<u32> pi((size_t)M + 1, 0);
    for (u64 i = 1; i <= M; ++i) pi[(size_t)i] = pi[(size_t)i - 1] + (is_prime[(size_t)i] ? 1U : 0U);

    std::vector<u64> cnt(32, 0); // p(n,k)

    for (u64 m = 1; m <= M; ++m) {
        const u64 pm = (u64)primes[(size_t)m - 1];
        const u64 A = (m < M) ? ((u64)primes[(size_t)m] - pm) : (n - pm + 1); // |{u0<=n: pi(u0)=m}|
        const u64 composites = (A >= 1) ? (A - 1) : 0;

        u64 u = m;
        int c = 0;
        while (u >= 1) {
            if (!is_prime[(size_t)u]) ++c;
            if ((size_t)(c + 2) > cnt.size()) cnt.resize((size_t)(c + 2), 0);
            cnt[(size_t)c] += 1;               // u0 is the unique prime in this pi-interval
            cnt[(size_t)(c + 1)] += composites; // u0 composite shifts c by +1
            if (u == 1) break;
            u = pi[(size_t)u];
        }
    }

    u64 prod = 1;
    for (u64 x : cnt) {
        if (x == 0) continue;
        prod = (u64)((u128)prod * (x % MOD) % MOD);
    }
    return prod;
}

int main() {
    // Statement validations.
    if (P(10) != 648ULL) {
        std::cerr << "Validation failed: P(10)\n";
        return 1;
    }
    if (P(100) != 31038676032ULL % MOD) {
        // P(100) fits in 64-bit; compare modulo MOD.
        std::cerr << "Validation failed: P(100)\n";
        return 1;
    }

    std::cout << P(100000000ULL) << "\n";
    return 0;
}
