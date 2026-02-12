#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;

static std::vector<int> primes_up_to(int n) {
    if (n < 2) {
        return {};
    }

    std::vector<char> is_prime(n + 1, 1);
    is_prime[0] = 0;
    is_prime[1] = 0;

    const int r = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    for (int p = 2; p <= r; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        const i64 start = 1LL * p * p;
        for (i64 x = start; x <= n; x += p) {
            is_prime[static_cast<int>(x)] = 0;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<int>(n / std::log(n)) + 10);
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

static i64 count_hybrid(int base, int exponent) {
    const long double limit = static_cast<long double>(exponent) * std::log(static_cast<long double>(base));
    const int q_max = static_cast<int>(limit / std::log(2.0L)) + 16;

    const std::vector<int> primes = primes_up_to(q_max);
    const int m = static_cast<int>(primes.size());

    std::vector<long double> logs(m);
    for (int i = 0; i < m; ++i) {
        logs[i] = std::log(static_cast<long double>(primes[i]));
    }

    i64 ans = 0;
    int j = m - 1;

    for (int i = 0; i < m; ++i) {
        if (i >= j) {
            break;
        }

        const int p = primes[i];
        const long double lp = logs[i];

        while (i < j) {
            const int q = primes[j];
            const long double lhs = static_cast<long double>(q) * lp +
                                    static_cast<long double>(p) * logs[j];
            if (lhs <= limit + 1e-18L) {
                break;
            }
            --j;
        }

        if (j <= i) {
            break;
        }
        ans += static_cast<i64>(j - i);
    }

    return ans;
}

int main() {
    assert(count_hybrid(800, 1) == 2);
    assert(count_hybrid(800, 800) == 10'790);

    std::cout << count_hybrid(800'800, 800'800) << '\n';
    return 0;
}
