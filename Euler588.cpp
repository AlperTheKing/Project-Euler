#include <array>
#include <cstdint>
#include <iostream>

// Project Euler 588: Quintinomial Coefficients
//
// Let P(x) = 1 + x + x^2 + x^3 + x^4 over F2.
// Q(k) is the number of odd coefficients in P(x)^k, i.e. the number of nonzero coefficients
// in P(x)^k modulo 2.
//
// Mod 2, P(x)^k is the generating function for sums of the form:
//   m = sum_{b>=0} d_b 2^b,
// where for each bit position b:
//   if the b-th bit of k is 0 then d_b = 0,
//   if the b-th bit of k is 1 then d_b can be 0..4.
// Each coefficient parity is the parity of the number of such digit assignments producing m.
//
// We compute Q(k) by running a carry-DP automaton over the bits of m, counting how many m
// yield odd parity. Carries are in {0,1,2,3,4} since the digit maximum is 4.
// The DP state is a 5-bit mask encoding the parity vector dp[carry].
// For each output bit (m_b) and for each bit of k (0 or 1) we get a deterministic transition.
// Counting accepted outputs is then a simple DP over the 32 masks.
//
// Also, since Q(2n) = Q(n) (squaring over F2 doubles exponents), Q(10^t) = Q(5^t), but we
// compute Q(10^t) directly.

using u64 = std::uint64_t;
using u8 = std::uint8_t;
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

static int bit_length_u64(u64 x) {
    if (x == 0) return 0;
    return 64 - __builtin_clzll(x);
}

static u128 Q(u64 k) {
    // next[t][out_bit][mask] -> new_mask
    static bool inited = false;
    static std::array<std::array<std::array<u8, 32>, 2>, 2> next{};
    if (!inited) {
        for (int t = 0; t <= 1; ++t) {
            for (int outb = 0; outb <= 1; ++outb) {
                for (int mask = 0; mask < 32; ++mask) {
                    int nmask = 0;
                    for (int c = 0; c <= 4; ++c) {
                        if (((mask >> c) & 1) == 0) continue;
                        const int dmax = t ? 4 : 0;
                        for (int d = 0; d <= dmax; ++d) {
                            const int s = c + d;
                            if ((s & 1) != outb) continue;
                            const int cp = s >> 1;
                            nmask ^= (1 << cp);
                        }
                    }
                    next[t][outb][mask] = (u8)nmask;
                }
            }
        }
        inited = true;
    }

    const int L = bit_length_u64(k) + 3; // enough to flush carry; max sum is 4k

    std::array<u128, 32> cnt{};
    cnt.fill(0);
    cnt[1] = 1; // dp[0]=1

    for (int b = 0; b < L; ++b) {
        const int t = (int)((k >> b) & 1ULL);
        std::array<u128, 32> nxt{};
        nxt.fill(0);
        for (int mask = 0; mask < 32; ++mask) {
            const u128 ways = cnt[mask];
            if (ways == 0) continue;
            nxt[next[t][0][mask]] += ways;
            nxt[next[t][1][mask]] += ways;
        }
        cnt = nxt;
    }

    u128 ans = 0;
    for (int mask = 0; mask < 32; ++mask) {
        if (mask & 1) ans += cnt[mask];
    }
    return ans;
}

int main() {
    // Validation points from the statement.
    if (Q(3) != 7) {
        std::cerr << "Validation failed: Q(3)\n";
        return 1;
    }
    if (Q(10) != 17) {
        std::cerr << "Validation failed: Q(10)\n";
        return 1;
    }
    if (Q(100) != 35) {
        std::cerr << "Validation failed: Q(100)\n";
        return 1;
    }

    u128 sum = 0;
    u64 k = 1;
    for (int e = 1; e <= 18; ++e) {
        k *= 10;
        sum += Q(k);
    }

    print_u128(sum);
    std::cout << "\n";
    return 0;
}
