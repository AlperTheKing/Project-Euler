#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

// Project Euler 617: count MPS by enumerating n = a + a^(e^L) and the valid starting roots of a.

using u64 = std::uint64_t;

static bool pow_leq(u64 a, int e, u64 lim) {
    __int128 r = 1;
    for (int i = 0; i < e; ++i) {
        r *= a;
        if (r > lim) return false;
    }
    return true;
}

static u64 floor_root(u64 n, int e) {
    if (e <= 1) return n;
    u64 lo = 1;
    u64 hi = std::min<u64>(n, 1000000000ULL) + 1;  // for this problem n <= 1e9 in all calls
    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;
        if (pow_leq(mid, e, n))
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

static u64 max_a_for(u64 N, int p) {
    u64 lo = 1, hi = 1000000000ULL + 1;
    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;
        __int128 r = 1;
        for (int i = 0; i < p; ++i) {
            r *= mid;
            if (r > (__int128)N) break;
        }
        if (r + mid <= (__int128)N)
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

// Sum_{a=2..A} chain_len(a,e), where chain_len is the number of times one can take an integer e-th root (>1).
static u64 sum_chain(u64 A, int e) {
    if (A < 2) return 0;
    u64 total = A - 1;  // k=0 term
    int t = e;
    while (t <= 60) {
        u64 r = floor_root(A, t);
        if (r < 2) break;
        total += (r - 1);
        if (t > 60 / e) break;
        t *= e;
    }
    return total;
}

static u64 compute_D(u64 N) {
    const u64 Nm2 = N - 2;
    int p_max = 1;
    __int128 v = 2;
    while ((v << 1) <= (__int128)Nm2) {
        v <<= 1;
        ++p_max;
    }

    std::vector<u64> A((std::size_t)p_max + 1, 0);
    for (int p = 2; p <= p_max; ++p) A[p] = max_a_for(N, p);

    u64 ans = 0;
    // L=1, e=p
    for (int p = 2; p <= p_max; ++p) ans += sum_chain(A[p], p);

    // L>=2, p=e^L
    for (int e = 2; e <= p_max; ++e) {
        int p = e * e;
        int L = 2;
        while (p <= p_max) {
            u64 amax = A[p];
            if (amax >= 2) ans += (u64)(L - 1) * (amax - 1);
            ans += sum_chain(amax, e);
            if (p > p_max / e) break;
            p *= e;
            ++L;
        }
    }
    return ans;
}

int main() {
#ifndef NDEBUG
    assert(compute_D(10) == 2);
    assert(compute_D(100) == 21);
    assert(compute_D(1000) == 69);
    assert(compute_D(1000000) == 1303);
    assert(compute_D(1000000000000ULL) == 1014800);
#endif

    std::cout << compute_D(1000000000000000000ULL) << "\n";
    return 0;
}

