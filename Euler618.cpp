#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

// Project Euler 618: DP on the generating function prod_p 1/(1 - p*x^p) up to max Fibonacci index.

using u64 = std::uint64_t;

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime((std::size_t)n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; (long long)p * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (int j = p * p; j <= n; j += p) is_prime[j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

int main() {
    static constexpr u64 MOD = 1000000000ULL;

    int F[25];
    F[1] = 1;
    F[2] = 1;
    for (int i = 3; i <= 24; ++i) F[i] = F[i - 1] + F[i - 2];
    const int K = F[24];

    std::vector<int> primes = sieve_primes(K);
    std::vector<u64> dp((std::size_t)K + 1, 0);
    dp[0] = 1;
    for (int p : primes) {
        for (int k = p; k <= K; ++k) {
            dp[k] += (u64)p * dp[k - p] % MOD;
            if (dp[k] >= MOD) dp[k] -= MOD;
        }
    }

#ifndef NDEBUG
    assert(dp[1] == 0);
    assert(dp[2] == 2);
    assert(dp[3] == 3);
    assert(dp[5] == 11);
    assert(dp[8] == 49);
#endif

    u64 ans = 0;
    for (int i = 2; i <= 24; ++i) {
        ans += dp[F[i]];
        ans %= MOD;
    }
    std::cout << std::setw(9) << std::setfill('0') << ans << "\n";
    return 0;
}

