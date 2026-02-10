#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;

static constexpr int MOD = 1'000'000'007;

static std::vector<int> primes_upto(int n) {
    std::vector<bool> is_comp(n / 2 + 1, false);
    std::vector<int> primes;
    primes.reserve(6'000'000);
    primes.push_back(2);
    for (int i = 3; (i64)i * i <= n; i += 2) {
        if (is_comp[i / 2]) continue;
        for (int j = i * i; j <= n; j += i * 2) is_comp[j / 2] = true;
    }
    for (int i = 3; i <= n; i += 2) {
        if (!is_comp[i / 2]) primes.push_back(i);
    }
    return primes;
}

static inline i64 isqrt_floor(i64 x) {
    i64 r = (i64)std::sqrt((long double)x);
    while ((r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

static std::array<u64, 10> compute_A(u64 N, const std::vector<uint8_t> &omega,
                                     const std::vector<uint8_t> &sqfree) {
    const int limit = (int)isqrt_floor((i64)N);
    std::array<u64, 10> A{};
    int i = 1;
    while (i <= limit) {
        const u64 q = N / (u64)i / (u64)i;
        int r = (int)isqrt_floor((i64)(N / q));
        if (r > limit) r = limit;
        while (r + 1 <= limit && N / (u64)(r + 1) / (u64)(r + 1) == q) ++r;
        while (r > limit || N / (u64)r / (u64)r != q) --r;

        for (int d = i; d <= r; ++d) {
            if (!sqfree[d]) continue;
            const int w = omega[d];
            if (w < 10) A[w] += q;
        }
        i = r + 1;
    }
    return A;
}

static std::array<u64, 10> compute_C(u64 N, const std::vector<uint8_t> &omega,
                                     const std::vector<uint8_t> &sqfree) {
    std::array<std::array<i64, 10>, 10> binom{};
    for (int n = 0; n < 10; ++n) {
        binom[n][0] = binom[n][n] = 1;
        for (int k = 1; k < n; ++k) binom[n][k] = binom[n - 1][k - 1] + binom[n - 1][k];
    }

    const auto A = compute_A(N, omega, sqfree);
    std::array<u64, 10> C{};
    for (int k = 0; k < 10; ++k) {
        __int128 s = 0;
        for (int r = k; r < 10; ++r) {
            const __int128 term = (__int128)A[r] * binom[r][k] * (((r - k) & 1) ? -1 : 1);
            s += term;
        }
        assert(s >= 0 && s <= (__int128)N);
        C[k] = (u64)s;
    }
    return C;
}

int main() {
    static constexpr int LIM = 100'000'000;
    const auto primes = primes_upto(LIM);

    std::vector<uint8_t> omega(LIM + 1, 0);
    std::vector<uint8_t> sqfree(LIM + 1, 1);
    sqfree[0] = 0;

    for (int p : primes) {
        for (int m = p; m <= LIM; m += p) ++omega[m];
        const i64 p2 = (i64)p * p;
        if (p2 <= LIM) {
            for (int m = (int)p2; m <= LIM; m += (int)p2) sqfree[m] = 0;
        }
    }

    {
        const auto C10 = compute_C(10ULL, omega, sqfree);
        assert(C10[0] == 7 && C10[1] == 3);
        for (int k = 2; k < 10; ++k) assert(C10[k] == 0);
    }
    {
        const auto C = compute_C(100ULL, omega, sqfree);
        assert(C[0] == 61 && C[1] == 36 && C[2] == 3);
    }
    {
        const auto C = compute_C(1'000ULL, omega, sqfree);
        assert(C[0] == 608 && C[1] == 343 && C[2] == 48 && C[3] == 1);
    }
    {
        const auto C = compute_C(10'000ULL, omega, sqfree);
        assert(C[0] == 6083 && C[1] == 3363 && C[2] == 533 && C[3] == 21);
    }
    {
        const auto C = compute_C(100'000ULL, omega, sqfree);
        assert(C[0] == 60794 && C[1] == 33562 && C[2] == 5345 && C[3] == 297 && C[4] == 2);
    }
    {
        const auto C = compute_C(1'000'000ULL, omega, sqfree);
        assert(C[0] == 607926 && C[1] == 335438 && C[2] == 53358 && C[3] == 3218 && C[4] == 60);
    }
    {
        const auto C = compute_C(10'000'000ULL, omega, sqfree);
        assert(C[0] == 6079291 && C[1] == 3353956 && C[2] == 533140 && C[3] == 32777 && C[4] == 834 &&
               C[5] == 2);
    }
    {
        const auto C = compute_C(100'000'000ULL, omega, sqfree);
        assert(C[0] == 60792694 && C[1] == 33539196 && C[2] == 5329747 && C[3] == 329028 && C[4] == 9257 &&
               C[5] == 78);
    }

    const u64 N = 10'000'000'000'000'000ULL;
    const auto C = compute_C(N, omega, sqfree);
    {
        __int128 s = 0;
        for (u64 x : C) s += x;
        assert(s == (__int128)N);
    }

    i64 ans = 1;
    for (u64 x : C) {
        if (x == 0) continue;
        ans = (i64)((__int128)ans * (x % MOD) % MOD);
    }
    std::cout << ans << "\n";
    return 0;
}

