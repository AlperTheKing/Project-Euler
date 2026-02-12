#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

static long double F_formula(long double n) {
    return (7.0L * n + 15.0L) / (18.0L * (n + 1.0L));
}

static std::vector<int> primes_upto(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int p = 2; static_cast<long long>(p) * p <= n; ++p) {
        if (is_prime[p]) {
            for (int q = p * p; q <= n; q += p) {
                is_prime[q] = false;
            }
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

int main() {
    assert(std::fabsl(F_formula(3.0L) - 0.5L) < 1e-18L);

    const int limit = 1'000'000;
    const auto primes = primes_upto(limit - 1);

    long double sum = 0.0L;
    std::uint64_t cnt = 0;

    for (int p : primes) {
        if (p == 2) {
            continue;
        }
        sum += F_formula(static_cast<long double>(p));
        ++cnt;
    }

    const long double ans = sum / static_cast<long double>(cnt);
    std::cout << std::fixed << std::setprecision(10) << ans << '\n';
    return 0;
}
