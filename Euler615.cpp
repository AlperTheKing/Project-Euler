#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 TARGET_COUNT = 1'000'000ULL;
constexpr u64 MOD = 123'454'321ULL;
constexpr int PRIME_LIMIT = 1'000'000;

u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    while (e > 0) {
        if (e & 1ULL) r = static_cast<u64>((static_cast<u128>(r) * a) % mod);
        a = static_cast<u64>((static_cast<u128>(a) * a) % mod);
        e >>= 1ULL;
    }
    return r;
}

std::vector<int> sieve_primes(int n) {
    std::vector<char> is_prime(static_cast<size_t>(n + 1), 1);
    is_prime[0] = 0;
    is_prime[1] = 0;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = 0;
    }
    std::vector<int> primes;
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

struct Enumerator {
    const std::vector<int>& primes;
    const std::vector<long double>& log_primes;
    u64 need_factors;
    u64 limit_n;
    long double log_limit_n;
    std::vector<u64> values;

    bool dfs(size_t i, u64 used, u64 n) {
        if (i >= primes.size()) return false;
        if (used < need_factors) {
            const long double lb = std::log(static_cast<long double>(n)) +
                                   static_cast<long double>(need_factors - used) * log_primes[i];
            if (lb > log_limit_n + 1e-18L) return false;
        }
        if (used >= need_factors) values.push_back(n);

        size_t j = i;
        while (j < primes.size()) {
            const u64 p = static_cast<u64>(primes[j]);
            if (static_cast<u128>(n) * p > limit_n) break;
            if (!dfs(j, used + 1, n * p)) break;
            ++j;
        }
        return true;
    }
};

u64 brute_sample() {
    auto omega = [](u64 n) {
        u64 cnt = 0;
        for (u64 p = 2; p * p <= n; ++p) {
            while (n % p == 0) {
                n /= p;
                ++cnt;
            }
        }
        if (n > 1) ++cnt;
        return cnt;
    };
    std::vector<u64> a;
    for (u64 n = 2; n <= 200; ++n) {
        if (omega(n) >= 5) a.push_back(n);
    }
    std::sort(a.begin(), a.end());
    return a[4];
}

u64 solve() {
    const std::vector<int> primes = sieve_primes(PRIME_LIMIT);
    std::vector<long double> log_primes(primes.size());
    for (size_t i = 0; i < primes.size(); ++i) {
        log_primes[i] = std::log(static_cast<long double>(primes[i]));
    }

    u64 c = 1;
    u64 n_limit = 3;
    std::vector<u64> values;

    while (true) {
        Enumerator e{primes, log_primes, c, n_limit, std::log(static_cast<long double>(n_limit)), {}};
        e.values.reserve(static_cast<size_t>(TARGET_COUNT + 1000));
        e.dfs(0, 0, 1);
        values.swap(e.values);
        if (values.size() >= TARGET_COUNT) break;
        ++c;
        n_limit *= 3;
    }

    std::sort(values.begin(), values.end());
    u64 x = values[static_cast<size_t>(TARGET_COUNT - 1)];
    for (u64 t = c; t < TARGET_COUNT; ++t) {
        x = static_cast<u64>((static_cast<u128>(x) * 2ULL) % MOD);
    }
    return x;
}

int main() {
    assert(brute_sample() == 80ULL);
    std::cout << solve() << '\n';
    return 0;
}
