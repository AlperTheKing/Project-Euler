#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

// Project Euler 616: creative n <= 1e12 are exactly the perfect powers except p^q (p,q primes) and 16.

using i64 = long long;
using u64 = std::uint64_t;

static bool is_prime_int(int x) {
    if (x < 2) return false;
    if (x % 2 == 0) return x == 2;
    for (int d = 3; (i64)d * d <= x; d += 2)
        if (x % d == 0) return false;
    return true;
}

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime((std::size_t)n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int p = 2; (i64)p * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (int j = p * p; j <= n; j += p) is_prime[j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

static bool pow_leq(i64 a, int e, i64 limit) {
    __int128 r = 1;
    for (int i = 0; i < e; ++i) {
        r *= a;
        if (r > limit) return false;
    }
    return true;
}

static i64 pow_exact(i64 a, int e) {
    __int128 r = 1;
    for (int i = 0; i < e; ++i) r *= a;
    assert(r <= std::numeric_limits<i64>::max());
    return (i64)r;
}

static i64 floor_root(i64 n, int e) {
    i64 lo = 1, hi = 1000000 + 1;  // since n <= 1e12 and e >= 2
    while (lo + 1 < hi) {
        i64 mid = lo + (hi - lo) / 2;
        if (pow_leq(mid, e, n))
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

static u64 pow_u64(u64 a, u64 e) {
    __int128 r = 1;
    for (u64 i = 0; i < e; ++i) {
        r *= a;
        assert(r <= (__int128)std::numeric_limits<u64>::max());
    }
    return (u64)r;
}

#ifndef NDEBUG
static void ms_erase_one(std::vector<u64> &L, u64 x) {
    auto it = std::find(L.begin(), L.end(), x);
    assert(it != L.end());
    L.erase(it);
}

static void op_split(std::vector<u64> &L, u64 c, u64 a, u64 b) {
    assert(a > 1 && b > 1);
    assert(pow_u64(a, b) == c);
    ms_erase_one(L, c);
    L.push_back(a);
    L.push_back(b);
}

static void op_merge(std::vector<u64> &L, u64 a, u64 b) {
    assert(a > 1 && b > 1);
    ms_erase_one(L, a);
    ms_erase_one(L, b);
    L.push_back(pow_u64(a, b));
}

static bool contains(const std::vector<u64> &L, u64 x) {
    return std::find(L.begin(), L.end(), x) != L.end();
}

static void validate() {
    // A few short constructive checks that "included" values can reach 64.
    {
        std::vector<u64> L{36};
        op_split(L, 36, 6, 2);
        op_merge(L, 2, 6);
        assert(contains(L, 64));
    }
    {
        std::vector<u64> L{81};
        op_split(L, 81, 3, 4);
        op_merge(L, 4, 3);
        assert(contains(L, 64));
    }
    {
        std::vector<u64> L{100};
        op_split(L, 100, 10, 2);
        op_merge(L, 2, 10);  // 2^10 = 1024
        op_split(L, 1024, 32, 2);
        op_merge(L, 2, 32);  // 2^32
        op_split(L, 4294967296ULL, 16, 8);
        op_split(L, 8, 2, 3);
        op_split(L, 16, 4, 2);
        op_merge(L, 4, 3);
        assert(contains(L, 64));
    }
    {
        std::vector<u64> L{216};
        op_split(L, 216, 6, 3);
        op_merge(L, 3, 6);  // 3^6
        op_split(L, 729, 27, 2);
        op_merge(L, 2, 27);  // 2^27
        op_split(L, 134217728ULL, 8, 9);
        op_split(L, 8, 2, 3);
        op_split(L, 9, 3, 2);
        op_merge(L, 2, 2);
        op_merge(L, 4, 3);
        assert(contains(L, 64));
    }

    // A tiny sanity check that 16 cannot introduce 3 (all operations stay within powers of 2 here).
    {
        std::vector<u64> L{16};
        op_split(L, 16, 2, 4);
        op_split(L, 4, 2, 2);
        // Any merge of powers of 2 stays a power of 2, and remaining elements are powers of 2 too.
        assert(!contains(L, 3));
    }
}
#endif

int main() {
    constexpr i64 LIM = 1000000000000LL;

#ifndef NDEBUG
    validate();
#endif

    int max_e = 1;
    while (pow_leq(2, max_e + 1, LIM)) ++max_e;

    std::vector<i64> perfect;
    perfect.reserve(1100000);
    for (int e = 2; e <= max_e; ++e) {
        i64 max_a = floor_root(LIM, e);
        for (i64 a = 2; a <= max_a; ++a) perfect.push_back(pow_exact(a, e));
    }
    std::sort(perfect.begin(), perfect.end());
    perfect.erase(std::unique(perfect.begin(), perfect.end()), perfect.end());

    u64 sum_perfect = 0;
    for (i64 v : perfect) sum_perfect += (u64)v;

    const int P_MAX = 1000000;
    std::vector<int> primes = sieve_primes(P_MAX);
    std::vector<int> exp_primes;
    for (int e = 2; e <= max_e; ++e)
        if (is_prime_int(e)) exp_primes.push_back(e);

    u64 sum_primeprime = 0;
    for (int p : primes) {
        for (int e : exp_primes) {
            if (!pow_leq(p, e, LIM)) break;
            sum_primeprime += (u64)pow_exact(p, e);
        }
    }

    u64 ans = sum_perfect - sum_primeprime - 16ULL;
    std::cout << ans << "\n";
    return 0;
}
