#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

using i64 = std::int64_t;
using i128 = __int128_t;

static std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }
    bool neg = false;
    if (x < 0) {
        neg = true;
        x = -x;
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    if (neg) {
        s.push_back('-');
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static i64 g_bruteforce(int n) {
    i64 sum = 0;
    for (int i = 1; i <= n; ++i) {
        const i64 term = std::gcd<i64>(n, 1LL * i * i);
        if (i & 1) {
            sum -= term;
        } else {
            sum += term;
        }
    }
    return sum;
}

static i128 G_bruteforce(int N) {
    i128 sum = 0;
    for (int n = 1; n <= N; ++n) {
        sum += g_bruteforce(n);
    }
    return sum;
}

static i128 solve(int N) {
    std::vector<int> spf(N + 1, 0);
    std::vector<int> primes;
    primes.reserve(N / 10);

    for (int i = 2; i <= N; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = 1LL * p * i;
            if (v > N) {
                break;
            }
            spf[static_cast<int>(v)] = p;
            if (p == spf[i]) {
                break;
            }
        }
    }

    std::vector<unsigned char> exp(N + 1, 0);
    std::vector<i64> c(N + 1, 0);
    c[1] = 1;

    for (int n = 2; n <= N; ++n) {
        const int p = spf[n];
        const int m = n / p;

        if (m % p == 0) {
            const int e = static_cast<int>(exp[m]) + 1;
            exp[n] = static_cast<unsigned char>(e);
            if (e % 2 == 0) {
                c[n] = c[m] * 1LL * p * p;
            } else {
                c[n] = c[m] * 1LL * p;
            }
        } else {
            exp[n] = 1;
            c[n] = c[m] * static_cast<i64>(p - 1);
        }
    }

    const i64 odd_count = (N + 1LL) / 2LL;
    i128 total = -static_cast<i128>(odd_count) * static_cast<i128>(odd_count);

    for (int d = 2; d <= N; d += 2) {
        const i64 m = N / d;
        const i128 tri = static_cast<i128>(m) * static_cast<i128>(m + 1) / 2;
        total += static_cast<i128>(c[d]) * tri;
    }

    return total;
}

int main() {
    assert(g_bruteforce(4) == 6);
    assert(g_bruteforce(1234) == 1233);
    assert(solve(1234) == static_cast<i128>(2'194'708));
    assert(solve(120) == G_bruteforce(120));

    std::cout << to_string_i128(solve(12'345'678)) << '\n';
    return 0;
}
