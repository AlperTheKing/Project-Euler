#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 clmul(u64 x, u64 y) {
    u64 r = 0;
    while (y > 0) {
        if (y & 1ULL) r ^= x;
        x <<= 1ULL;
        y >>= 1ULL;
    }
    return r;
}

static u64 value_f(u64 a, u64 b) {
    return clmul(a, a) ^ clmul(clmul(2, a), b) ^ clmul(b, b);
}

static u64 seed_bound(u64 m) {
    if (m == 0) return 1;
    unsigned bits = std::bit_width(m);
    return 1ULL << ((bits + 2) / 2);
}

static u64 count_fast(u64 n, u64 m) {
    if (m == 0) return 1;

    u64 bmax = seed_bound(m);
    u64 total = 0;

    for (u64 b = 0; b < bmax; ++b) {
        for (u64 a = 0; a <= b; ++a) {
            u64 k = value_f(a, b);
            if (k > m) continue;

            if (a == 0 && b == 0) {
                ++total;
                continue;
            }

            u64 prev = b ^ (a << 1ULL);
            if (prev <= a) continue;

            u128 aa = a;
            u128 bb = b;
            while (bb <= n) {
                ++total;
                u128 next = (bb << 1ULL) ^ aa;
                aa = bb;
                bb = next;
            }
        }
    }

    return total;
}

static u64 count_bruteforce(u64 n, u64 m) {
    u64 total = 0;
    for (u64 b = 0; b <= n; ++b) {
        for (u64 a = 0; a <= b; ++a) {
            if (value_f(a, b) <= m) ++total;
        }
    }
    return total;
}

int main() {
    assert(count_fast(1000, 100) == 398ULL);

    for (u64 n = 0; n <= 80; n += 10) {
        for (u64 m = 0; m <= 200; m += 25) {
            assert(count_fast(n, m) == count_bruteforce(n, m));
        }
    }

    std::cout << count_fast(100'000'000'000'000'000ULL, 1'000'000ULL) << '\n';
    return 0;
}
