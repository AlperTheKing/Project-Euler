#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_dec_float_100;

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = false;
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[i]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[j] = false;
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

static long double brute_f1_average(int n) {
    std::vector<int> spf(n + 1, 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[i] != 0) {
            continue;
        }
        spf[i] = i;
        if (1LL * i * i > n) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            if (spf[j] == 0) {
                spf[j] = i;
            }
        }
    }
    long double sum = 0.0L;
    for (int x = 2; x <= n; ++x) {
        int p = spf[x];
        int t = x;
        int a = 0;
        while (t % p == 0) {
            t /= p;
            ++a;
        }
        sum += static_cast<long double>(a - 1) / static_cast<long double>(p);
    }
    return sum / static_cast<long double>(n - 1);
}

int main() {
    using Real = cpp_dec_float_100;

    const int max_prime = 2'000'000;
    const auto primes = sieve_primes(max_prime);

    Real prefix = 1;
    Real f1 = 0;
    Real total = 0;

    for (int p_int : primes) {
        const Real p = Real(p_int);
        const Real pm1 = p - 1;

        f1 += prefix / (p * p * pm1);
        total += prefix / (p * pm1 * pm1);
        prefix *= pm1 / p;
    }

    const Real tail_bound = Real(1) / (Real(2) * Real(max_prime) * Real(max_prime));
    assert(tail_bound < Real("5e-13"));

    const Real stated_f1 = Real("0.282419756159");
    assert(abs(f1 - stated_f1) < Real("1e-9"));

    const long double brute = brute_f1_average(200000);
    const long double theo = f1.convert_to<long double>();
    assert(std::fabsl(brute - theo) < 2e-4L);

    std::cout << std::fixed << std::setprecision(12) << total << '\n';
    return 0;
}
