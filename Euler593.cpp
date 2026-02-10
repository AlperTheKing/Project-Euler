#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

// Project Euler 593: Fleeting Medians
//
// S(k)  = p_k^k mod 10007, where p_k is the k-th prime.
// S2(k) = S(k) + S(floor(k/10000) + 1).
//
// Values are small: S(k) in [0,10006], so S2(k) in [0,20012].
// For sliding-window medians we maintain a Fenwick tree of frequencies and query order
// statistics in O(log V) time.
//
// Sum of medians can be half-integer; we accumulate twice the answer to stay integral.

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }
    char buf[64];
    int n = 0;
    while (x > 0) {
        const u64 digit = (u64)(x % 10);
        buf[n++] = (char)('0' + digit);
        x /= 10;
    }
    while (n--) std::cout << buf[n];
}

static u64 nth_prime_upper_bound(u64 n) {
    if (n < 6) return 15;
    const long double nn = (long double)n;
    const long double bound = nn * (logl(nn) + logl(logl(nn)));
    return (u64)(bound + 16.0L);
}

struct Fenwick {
    int n = 0;
    std::vector<int> bit;

    explicit Fenwick(int n_) : n(n_), bit(n_ + 1, 0) {}

    void add(int idx0, int delta) {
        // idx0 is 0-based value; Fenwick is 1-based.
        int i = idx0 + 1;
        for (; i <= n; i += i & -i) bit[i] += delta;
    }

    int kth(int k) const {
        // Smallest value v (0-based) such that prefix sum >= k, assuming 1 <= k <= total.
        int idx = 0;
        int step = 1;
        while (step << 1 <= n) step <<= 1;
        for (; step; step >>= 1) {
            const int nxt = idx + step;
            if (nxt <= n && bit[nxt] < k) {
                idx = nxt;
                k -= bit[nxt];
            }
        }
        return idx; // idx is 0-based because we return (idx+1)-1
    }
};

static int powmod_int(int a, int e, int mod) {
    int r = 1;
    int x = a % mod;
    while (e > 0) {
        if (e & 1) r = (int)((1LL * r * x) % mod);
        x = (int)((1LL * x * x) % mod);
        e >>= 1;
    }
    return r;
}

struct OddCompositeBits {
    // Composite flags for odd numbers up to limit.
    // Index is x>>1 (so 3 -> 1, 5 -> 2, ...). Bit=1 means composite.
    u64 limit;
    std::vector<u64> bits;

    explicit OddCompositeBits(u64 limit_) : limit(limit_) {
        const u64 sz_bits = (limit >> 1) + 1;
        bits.assign((sz_bits + 63) >> 6, 0ULL);
    }

    inline bool get(u64 half) const {
        return (bits[half >> 6] >> (half & 63ULL)) & 1ULL;
    }

    inline void set(u64 half) {
        bits[half >> 6] |= 1ULL << (half & 63ULL);
    }
};

static void sieve_odd(OddCompositeBits& comp) {
    const u64 limit = comp.limit;
    const u64 r = (u64)std::sqrt((long double)limit);
    for (u64 p = 3; p <= r; p += 2) {
        if (comp.get(p >> 1)) continue;
        const u64 step = p << 1;
        for (u64 x = p * p; x <= limit; x += step) comp.set(x >> 1);
    }
}

static std::vector<int> generate_S2(u64 n) {
    constexpr int MOD = 10007;
    constexpr int PHI = MOD - 1;

    const u64 limit = nth_prime_upper_bound(n);
    OddCompositeBits comp(limit);
    sieve_odd(comp);

    std::vector<int> Sfirst(1002, 0);
    std::vector<int> S2(n + 1, 0);

    u64 k = 0;

    auto handle_prime = [&](u64 p) {
        ++k;
        const int base = (int)(p % MOD);
        int s = 0;
        if (base != 0) {
            const int e = (int)(k % PHI);
            s = powmod_int(base, e, MOD);
        }

        if (k <= 1001) Sfirst[(size_t)k] = s;
        const int idx = (int)(k / 10000 + 1);
        const int s2 = s + Sfirst[(size_t)idx];
        S2[(size_t)k] = s2;
    };

    handle_prime(2);
    for (u64 x = 3; x <= limit && k < n; x += 2) {
        if (!comp.get(x >> 1)) handle_prime(x);
    }

    assert(k == n);
    return S2;
}

static u128 median2_from_sorted(const std::vector<int>& v) {
    const size_t len = v.size();
    if (len & 1) return (u128)2 * (u128)v[len / 2];
    return (u128)v[len / 2 - 1] + (u128)v[len / 2];
}

static u128 median2_range(const std::vector<int>& S2, u64 l, u64 r) {
    std::vector<int> tmp;
    tmp.reserve((size_t)(r - l + 1));
    for (u64 i = l; i <= r; ++i) tmp.push_back(S2[(size_t)i]);
    std::sort(tmp.begin(), tmp.end());
    return median2_from_sorted(tmp);
}

static u128 F2(u64 n, int K) {
    constexpr int MOD = 10007;
    constexpr int PHI = MOD - 1;
    constexpr int VMAX = 20013; // values 0..20012

    const u64 limit = nth_prime_upper_bound(n);
    OddCompositeBits comp(limit);
    sieve_odd(comp);

    std::vector<int> Sfirst(1002, 0);
    Fenwick fw(VMAX);
    std::vector<std::uint16_t> ring((size_t)K);
    int ptr = 0;

    u128 sum2 = 0;

    auto add_median = [&]() {
        if ((K & 1) == 1) {
            const int v = fw.kth((K + 1) / 2);
            sum2 += (u128)2 * (u128)v;
        } else {
            const int a = fw.kth(K / 2);
            const int b = fw.kth(K / 2 + 1);
            sum2 += (u128)a + (u128)b;
        }
    };

    u64 k = 0;

    auto handle_prime = [&](u64 p) {
        ++k;
        const int base = (int)(p % MOD);
        int s = 0;
        if (base != 0) {
            const int e = (int)(k % PHI);
            s = powmod_int(base, e, MOD);
        }
        if (k <= 1001) Sfirst[(size_t)k] = s;
        const int idx = (int)(k / 10000 + 1);
        const int v = s + Sfirst[(size_t)idx];

        if ((int)k <= K) {
            ring[(size_t)ptr++] = (std::uint16_t)v;
            fw.add(v, +1);
            if ((int)k == K) {
                ptr = 0;
                add_median();
            }
        } else {
            const int out = (int)ring[(size_t)ptr];
            fw.add(out, -1);
            ring[(size_t)ptr] = (std::uint16_t)v;
            fw.add(v, +1);
            ++ptr;
            if (ptr == K) ptr = 0;
            add_median();
        }
    };

    handle_prime(2);
    for (u64 x = 3; x <= limit && k < n; x += 2) {
        if (!comp.get(x >> 1)) handle_prime(x);
    }
    assert(k == n);

    return sum2;
}

int main() {
    // Median validations from the statement.
    {
        const auto S2 = generate_S2(1000);
        if (median2_range(S2, 1, 10) != 4043) {
            std::cerr << "Validation failed: M(1,10)\n";
            return 1;
        }
        if (median2_range(S2, 100, 1000) != 9430) {
            std::cerr << "Validation failed: M(100,1000)\n";
            return 1;
        }
    }

    // F validations from the statement.
    if (F2(100, 10) != 927257) {
        std::cerr << "Validation failed: F(100,10)\n";
        return 1;
    }
    if (F2(100000, 10000) != 1350696415ULL) {
        std::cerr << "Validation failed: F(1e5,1e4)\n";
        return 1;
    }

    const u128 ans2 = F2(10000000ULL, 100000);
    print_u128(ans2 / 2);
    std::cout << ((ans2 & 1) ? ".5\n" : ".0\n");
    return 0;
}
