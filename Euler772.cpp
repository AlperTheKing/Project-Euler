#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

constexpr i64 MOD = 1'000'000'007LL;

std::vector<int> odd_sieve_primes(const int n) {
    std::vector<int> primes;
    if (n < 2) {
        return primes;
    }
    primes.push_back(2);
    if (n < 3) {
        return primes;
    }

    const int m = (n - 1) / 2;  // odd numbers 3,5,...,2m+1
    std::vector<bool> composite(static_cast<std::size_t>(m), false);

    for (int i = 0; (2LL * i + 3) * (2LL * i + 3) <= n; ++i) {
        if (composite[static_cast<std::size_t>(i)]) {
            continue;
        }
        const int p = 2 * i + 3;
        i64 j = (static_cast<i64>(p) * p - 3) / 2;
        while (j < m) {
            composite[static_cast<std::size_t>(j)] = true;
            j += p;
        }
    }

    for (int i = 0; i < m; ++i) {
        if (!composite[static_cast<std::size_t>(i)]) {
            primes.push_back(2 * i + 3);
        }
    }
    return primes;
}

i64 lcm_upto_mod(const int k) {
    const std::vector<int> primes = odd_sieve_primes(k);
    i64 ans = 1;
    for (const int p : primes) {
        u64 pw = static_cast<u64>(p);
        while (pw <= static_cast<u64>(k) / static_cast<u64>(p)) {
            pw *= static_cast<u64>(p);
        }
        ans = static_cast<i64>((static_cast<__int128>(ans) * (pw % MOD)) % MOD);
    }
    return ans;
}

i64 f_mod(const int k) {
    const i64 l = lcm_upto_mod(k);
    return (2LL * l) % MOD;
}

}  // namespace

int main() {
    assert(f_mod(3) == 12);
    assert(f_mod(30) == 179'092'994LL);

    std::cout << f_mod(100'000'000) << '\n';
    return 0;
}
