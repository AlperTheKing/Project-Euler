#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

// Project Euler 615: generate the smallest products of {2, p/2 for odd primes p} and scale by 2^K.

using i64 = long long;
using u64 = std::uint64_t;

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e) {
        if (e & 1ULL) r = (u64)((__int128)r * a % mod);
        a = (u64)((__int128)a * a % mod);
        e >>= 1ULL;
    }
    return r;
}

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime((std::size_t)n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; (i64)p * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (int j = p * p; j <= n; j += p) is_prime[j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

struct Gen {
    long double logv;
    u64 mul1;
    u64 mul2;
};

struct Ratio {
    long double logv;
    u64 mod1;
    u64 mod2;
};

static Ratio millionth_ratio(const std::vector<int>& odd_primes, u64 mod1, u64 mod2) {
    const u64 inv2_1 = (mod1 + 1) / 2; // mod is odd
    const u64 inv2_2 = (mod2 + 1) / 2;
    const long double ln2 = std::log((long double)2.0);

    std::vector<Gen> gens;
    gens.reserve(odd_primes.size() + 1);
    gens.push_back(Gen{ln2, 2 % mod1, 2 % mod2}); // multiply by 2
    for (int p : odd_primes) {
        const long double lg = std::log((long double)p) - ln2;
        const u64 mm1 = (u64)((__int128)(p % (int)mod1) * inv2_1 % mod1);
        const u64 mm2 = (u64)((__int128)(p % (int)mod2) * inv2_2 % mod2);
        gens.push_back(Gen{lg, mm1, mm2}); // multiply by p/2
    }
    std::sort(gens.begin(), gens.end(), [](const Gen& a, const Gen& b) { return a.logv < b.logv; });

    const int M = (int)gens.size();
    const int K = 999'999; // 0-based index: h[0]=1 is the first number

    std::vector<long double> hlog((std::size_t)K + 1, 0.0L);
    std::vector<u64> h1((std::size_t)K + 1, 1ULL);
    std::vector<u64> h2((std::size_t)K + 1, 1ULL);
    std::vector<int> idx((std::size_t)M, 0);

    const long double eps = 1e-18L;
    for (int i = 1; i <= K; ++i) {
        long double best_log = std::numeric_limits<long double>::infinity();
        int best_j = -1;
        for (int j = 0; j < M; ++j) {
            const long double cand = hlog[idx[j]] + gens[j].logv;
            if (cand < best_log) {
                best_log = cand;
                best_j = j;
            }
        }
        assert(best_j >= 0);
        const u64 best1 = (u64)((__int128)h1[idx[best_j]] * gens[best_j].mul1 % mod1);
        const u64 best2 = (u64)((__int128)h2[idx[best_j]] * gens[best_j].mul2 % mod2);

        hlog[i] = best_log;
        h1[i] = best1;
        h2[i] = best2;

        for (int j = 0; j < M; ++j) {
            const long double cand = hlog[idx[j]] + gens[j].logv;
            if (fabsl(cand - best_log) > eps) continue;
            const u64 c1 = (u64)((__int128)h1[idx[j]] * gens[j].mul1 % mod1);
            const u64 c2 = (u64)((__int128)h2[idx[j]] * gens[j].mul2 % mod2);
            if (c1 == best1 && c2 == best2) ++idx[j];
        }
    }
    return Ratio{hlog[K], h1[K], h2[K]};
}

static void validate_small() {
    // K=5, the 5th number with Ω(n) >= 5 is 80.
    const u64 mod = 1'000'000'007ULL;
    const u64 inv2 = (mod + 1) / 2;
    const std::vector<Gen> gens = {
        {std::log(1.5L), (3ULL * inv2) % mod, (3ULL * inv2) % mod},
        {std::log(2.0L), 2ULL % mod, 2ULL % mod},
        {std::log(2.5L), (5ULL * inv2) % mod, (5ULL * inv2) % mod},
    };

    const int M = (int)gens.size();
    const int target = 4; // 0-based: 5th element is h[4]
    std::vector<long double> hlog((std::size_t)target + 1, 0.0L);
    std::vector<u64> hmod((std::size_t)target + 1, 1ULL);
    std::vector<int> idx((std::size_t)M, 0);

    const long double eps = 1e-18L;
    for (int i = 1; i <= target; ++i) {
        long double best_log = std::numeric_limits<long double>::infinity();
        int best_j = -1;
        for (int j = 0; j < M; ++j) {
            const long double cand = hlog[idx[j]] + gens[j].logv;
            if (cand < best_log) {
                best_log = cand;
                best_j = j;
            }
        }
        const u64 best_mod = (u64)((__int128)hmod[idx[best_j]] * gens[best_j].mul1 % mod);
        hlog[i] = best_log;
        hmod[i] = best_mod;
        for (int j = 0; j < M; ++j) {
            const long double cand = hlog[idx[j]] + gens[j].logv;
            if (fabsl(cand - best_log) > eps) continue;
            const u64 cand_mod = (u64)((__int128)hmod[idx[j]] * gens[j].mul1 % mod);
            if (cand_mod == best_mod) ++idx[j];
        }
    }

    const u64 n = (u64)((__int128)mod_pow(2, 5, mod) * hmod[target] % mod);
    assert(n == 80ULL);
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    validate_small();

    constexpr u64 K = 1'000'000ULL;
    constexpr u64 MODANS = 123'454'321ULL;
    constexpr u64 MODAUX = 1'000'000'007ULL;

    int limit = 200;
    Ratio r{0.0L, 0ULL, 0ULL};
    for (;;) {
        const auto primes = sieve_primes(limit);
        std::vector<int> odd_primes;
        for (int p : primes) {
            if (p >= 3) odd_primes.push_back(p);
        }
        r = millionth_ratio(odd_primes, MODANS, MODAUX);

        // Find the next odd prime after `limit`.
        int next_p = limit + 1;
        for (;;) {
            bool ok = true;
            for (int q = 2; (i64)q * q <= next_p; ++q) {
                if (next_p % q == 0) {
                    ok = false;
                    break;
                }
            }
            if (ok) break;
            ++next_p;
        }
        const long double ln2 = std::log((long double)2.0);
        const long double next_log = std::log((long double)next_p) - ln2;
        if (r.logv + 1e-15L < next_log) break;
        limit *= 2;
    }

    const u64 ratio_mod = r.mod1;
    const u64 base = mod_pow(2, K, MODANS);
    const u64 answer = (u64)((__int128)base * ratio_mod % MODANS);
    std::cout << answer << "\n";
    return 0;
}
