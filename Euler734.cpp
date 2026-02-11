#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr int kMod = 1'000'000'007;

int mod_pow(int base, int exp) {
    i64 b = base;
    i64 result = 1;
    while (exp > 0) {
        if ((exp & 1) != 0) {
            result = (result * b) % kMod;
        }
        b = (b * b) % kMod;
        exp >>= 1;
    }
    return static_cast<int>(result);
}

std::vector<int> sieve_primes(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; static_cast<i64>(p) * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int x = p * p; x <= n; x += p) {
            is_prime[static_cast<std::size_t>(x)] = false;
        }
    }
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    for (int x = 2; x <= n; ++x) {
        if (is_prime[static_cast<std::size_t>(x)]) {
            primes.push_back(x);
        }
    }
    return primes;
}

int T(const int n, const int k) {
    int bits = 0;
    while ((1 << bits) <= n) {
        ++bits;
    }
    const int size = 1 << bits;

    const std::vector<int> primes = sieve_primes(n);

    std::vector<int> f(static_cast<std::size_t>(size), 0);
    for (const int p : primes) {
        f[static_cast<std::size_t>(p)] = 1;
    }

    std::vector<int> g = f;
    for (int b = 0; b < bits; ++b) {
        const int bit = 1 << b;
        for (int mask = 0; mask < size; ++mask) {
            if ((mask & bit) == 0) {
                continue;
            }
            int v = g[static_cast<std::size_t>(mask)] + g[static_cast<std::size_t>(mask ^ bit)];
            if (v >= kMod) {
                v -= kMod;
            }
            g[static_cast<std::size_t>(mask)] = v;
        }
    }

    std::vector<int> exact(static_cast<std::size_t>(size), 0);
    for (int mask = 0; mask < size; ++mask) {
        exact[static_cast<std::size_t>(mask)] = mod_pow(g[static_cast<std::size_t>(mask)], k);
    }

    for (int b = 0; b < bits; ++b) {
        const int bit = 1 << b;
        for (int mask = 0; mask < size; ++mask) {
            if ((mask & bit) == 0) {
                continue;
            }
            int v = exact[static_cast<std::size_t>(mask)] - exact[static_cast<std::size_t>(mask ^ bit)];
            if (v < 0) {
                v += kMod;
            }
            exact[static_cast<std::size_t>(mask)] = v;
        }
    }

    int answer = 0;
    for (const int p : primes) {
        answer += exact[static_cast<std::size_t>(p)];
        if (answer >= kMod) {
            answer -= kMod;
        }
    }
    return answer;
}

}  // namespace

int main() {
    assert(T(100, 3) == 3355);
    assert(T(1000, 10) == 2'071'632);

    std::cout << T(1'000'000, 999'983) << '\n';
    return 0;
}
