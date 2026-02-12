#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;
using u32 = std::uint32_t;

static i64 floor_div(i64 a, i64 b) {
    if (a >= 0) {
        return a / b;
    }
    return -(((-a) + b - 1) / b);
}

static std::vector<u32> build_spf(int n) {
    std::vector<u32> spf(n + 1, 0);
    for (u32 i = 2; 1ULL * i * i <= static_cast<u32>(n); ++i) {
        if (spf[i] != 0) {
            continue;
        }
        for (u64 j = 1ULL * i * i; j <= static_cast<u64>(n); j += i) {
            if (spf[static_cast<std::size_t>(j)] == 0) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return spf;
}

static int coprime_prefix_count(int x, int t, const std::vector<u32>& spf) {
    if (t <= 0) {
        return 0;
    }

    int primes[10];
    int k = 0;
    int n = x;
    while (n > 1) {
        int p = static_cast<int>(spf[static_cast<std::size_t>(n)]);
        if (p == 0) {
            p = n;
        }
        primes[k++] = p;
        while (n % p == 0) {
            n /= p;
        }
    }

    int divs[1 << 10];
    int signs[1 << 10];
    int sz = 1;
    divs[0] = 1;
    signs[0] = 1;

    for (int i = 0; i < k; ++i) {
        const int p = primes[i];
        const int old = sz;
        for (int j = 0; j < old; ++j) {
            divs[sz] = divs[j] * p;
            signs[sz] = -signs[j];
            ++sz;
        }
    }

    int res = 0;
    for (int i = 0; i < sz; ++i) {
        res += signs[i] * (t / divs[i]);
    }
    return res;
}

static u64 solve_fast(i64 N) {
    if (N < 2) {
        return 0;
    }

    const i64 xmax1 = 3 * floor_div(N - 8, 10) + 1;
    const i64 xmax2 = 3 * floor_div(N - 11, 10) + 2;
    const int xmax = static_cast<int>(std::max<i64>(0, std::max(xmax1, xmax2)));

    const auto spf = build_spf(xmax);

    u64 sum = 0;
    for (int x = 1; x <= xmax; ++x) {
        if (x % 3 == 0) {
            continue;
        }

        const i64 t = (N + 1 - (x / 3) - 3LL * x) / 6;
        if (t < 1) {
            continue;
        }

        sum += static_cast<u64>(coprime_prefix_count(x, static_cast<int>(t), spf));
    }

    return 2ULL + 4ULL * sum;
}

static u64 solve_bruteforce(i64 N) {
    u64 ans = 0;
    const i64 mmax = N / 3 + 3;

    for (i64 m = 1; m <= mmax; ++m) {
        for (i64 n = (m & 1LL) ? 1LL : 2LL; n < 3 * m; n += 2) {
            const i64 v = (n - m) / 2;
            if (std::gcd(m, std::llabs(v)) != 1) {
                continue;
            }
            if ((m - v) % 3 == 0) {
                continue;
            }

            i64 bounces;
            if (v >= 0) {
                bounces = 3 * (m + v) - 1 + (m - v) / 3;
            } else {
                bounces = 3 * m - 1 + n / 3;
            }

            if (bounces <= N) {
                ans += 2;
            }
        }
    }

    return ans;
}

int main() {
    assert(solve_fast(10) == 6);
    assert(solve_fast(100) == 478);
    assert(solve_fast(1000) == 45790);
    assert(solve_fast(2000) == solve_bruteforce(2000));

    const u64 ans = solve_fast(1'000'000'000LL);
    std::cout << ans << '\n';
    return 0;
}
