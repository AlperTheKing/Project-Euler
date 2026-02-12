#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;

static inline int deg_u32(u32 x) {
    return x ? 31 - __builtin_clz(x) : -1;
}

static inline int deg_u64(u64 x) {
    return x ? 63 - __builtin_clzll(x) : -1;
}

static u32 poly_mod_u64(u64 a, u32 b) {
    const int db = deg_u32(b);
    while (a && deg_u64(a) >= db) {
        const int sh = deg_u64(a) - db;
        a ^= static_cast<u64>(b) << sh;
    }
    return static_cast<u32>(a);
}

static u32 poly_mod_u32(u32 a, u32 b) {
    const int db = deg_u32(b);
    while (a && deg_u32(a) >= db) {
        const int sh = deg_u32(a) - db;
        a ^= b << sh;
    }
    return a;
}

static u32 poly_gcd(u32 a, u32 b) {
    while (b) {
        const u32 r = poly_mod_u32(a, b);
        a = b;
        b = r;
    }
    return a;
}

static std::array<u32, 1u << 16> SQ16;

static void init_sq16() {
    for (u32 x = 0; x < (1u << 16); ++x) {
        u32 y = 0;
        for (u32 i = 0; i < 16; ++i) {
            if (x & (1u << i)) {
                y |= 1u << (2u * i);
            }
        }
        SQ16[x] = y;
    }
}

static inline u32 square_mod(u32 a, u32 mod_poly) {
    const u64 sq = static_cast<u64>(SQ16[a & 0xFFFFu]) | (static_cast<u64>(SQ16[a >> 16]) << 32);
    return poly_mod_u64(sq, mod_poly);
}

static std::vector<int> prime_divisors(int n) {
    std::vector<int> p;
    for (int d = 2; d * d <= n; ++d) {
        if (n % d == 0) {
            p.push_back(d);
            while (n % d == 0) {
                n /= d;
            }
        }
    }
    if (n > 1) {
        p.push_back(n);
    }
    return p;
}

static bool is_irreducible(u32 f, int d, const std::vector<int>& prime_factors) {
    if (d == 1) {
        return true;
    }

    u32 x = 2u;
    for (int i = 1; i <= d; ++i) {
        x = square_mod(x, f);
        for (int p : prime_factors) {
            if (i == d / p) {
                if (poly_gcd(f, x ^ 2u) != 1u) {
                    return false;
                }
            }
        }
    }

    return x == 2u;
}

static int mobius(int n) {
    int mu = 1;
    for (int p = 2; p * p <= n; ++p) {
        if (n % p == 0) {
            n /= p;
            if (n % p == 0) {
                return 0;
            }
            mu = -mu;
            while (n % p == 0) {
                n /= p;
            }
        }
    }
    if (n > 1) {
        mu = -mu;
    }
    return mu;
}

static std::vector<int> divisors(int n) {
    std::vector<int> d;
    for (int i = 1; i * i <= n; ++i) {
        if (n % i == 0) {
            d.push_back(i);
            if (i * i != n) {
                d.push_back(n / i);
            }
        }
    }
    std::sort(d.begin(), d.end());
    return d;
}

static u64 count_irreducibles_degree(int d) {
    long long s = 0;
    for (int k : divisors(d)) {
        s += static_cast<long long>(mobius(k)) * (1LL << (d / k));
    }
    return static_cast<u64>(s / d);
}

static u32 nth_xor_prime(u64 n) {
    if (n == 1) {
        return 2u;
    }
    if (n == 2) {
        return 3u;
    }

    u64 cumulative = 0;
    int target_degree = 0;
    u64 target_in_degree = 0;

    for (int d = 1;; ++d) {
        const u64 cnt = count_irreducibles_degree(d);
        if (cumulative + cnt >= n) {
            target_degree = d;
            target_in_degree = n - cumulative;
            break;
        }
        cumulative += cnt;
    }

    const std::vector<int> prime_factors = prime_divisors(target_degree);

    const u32 top = 1u << target_degree;
    const u32 mid_limit = 1u << (target_degree - 1);

    u64 seen = 0;
    for (u32 mid = 0; mid < mid_limit; ++mid) {
        const u32 f = top | (mid << 1) | 1u;
        if (is_irreducible(f, target_degree, prime_factors)) {
            ++seen;
            if (seen == target_in_degree) {
                return f;
            }
        }
    }

    return 0u;
}

int main() {
    init_sq16();

    assert(nth_xor_prime(10) == 41u);

    std::cout << nth_xor_prime(5'000'000ULL) << '\n';
    return 0;
}
