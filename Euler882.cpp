#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;

struct Dyad {
    i64 num;
    int exp;
};

static Dyad normalize(i128 num, int exp) {
    if (num == 0) return {0, 0};
    while (exp > 0 && (num & 1) == 0) {
        num >>= 1;
        --exp;
    }
    return {static_cast<i64>(num), exp};
}

static i128 floor_div_pow2(i128 x, int sh) {
    if (sh == 0) return x;
    if (x >= 0) return x >> sh;
    i128 a = -x;
    return -((a + ((static_cast<i128>(1) << sh) - 1)) >> sh);
}

static int cmp_dyad(const Dyad& a, const Dyad& b) {
    i128 na = a.num;
    i128 nb = b.num;
    if (a.exp < b.exp) na <<= (b.exp - a.exp);
    else if (b.exp < a.exp) nb <<= (a.exp - b.exp);
    if (na < nb) return -1;
    if (na > nb) return 1;
    return 0;
}

static i128 floor_dyad(const Dyad& a) {
    return floor_div_pow2(static_cast<i128>(a.num), a.exp);
}

static i128 ceil_dyad(const Dyad& a) {
    Dyad neg{-a.num, a.exp};
    return -floor_dyad(neg);
}

static i128 floor_scaled(const Dyad& a, int k) {
    i128 n = static_cast<i128>(a.num);
    if (k >= a.exp) return n << (k - a.exp);
    return floor_div_pow2(n, a.exp - k);
}

static i128 ceil_scaled(const Dyad& a, int k) {
    Dyad neg{-a.num, a.exp};
    return -floor_scaled(neg, k);
}

static Dyad simplest_between(const Dyad& L, bool hasL, const Dyad& R, bool hasR) {
    if (!hasL && !hasR) return {0, 0};
    if (hasL && !hasR) return {static_cast<i64>(floor_dyad(L) + 1), 0};
    if (!hasL && hasR) return {static_cast<i64>(ceil_dyad(R) - 1), 0};

    assert(cmp_dyad(L, R) < 0);

    for (int k = 0; k <= 80; ++k) {
        i128 lo = floor_scaled(L, k) + 1;
        i128 hi = ceil_scaled(R, k) - 1;
        if (lo > hi) continue;

        i128 m;
        if (lo <= 0 && hi >= 0) m = 0;
        else if (hi < 0) m = hi;
        else m = lo;

        return normalize(m, k);
    }

    assert(false);
    return {0, 0};
}

static int bit_length(i64 x) {
    if (x == 0) return 1;
    return 64 - __builtin_clzll(static_cast<unsigned long long>(x));
}

static i64 delete_bit(i64 x, int len, int idx_from_left) {
    int p = len - 1 - idx_from_left;
    i64 hi = x >> (p + 1);
    i64 lo = (p == 0) ? 0 : (x & ((1LL << p) - 1));
    return (hi << p) | lo;
}

static std::vector<Dyad> compute_values(int n) {
    std::vector<Dyad> v(static_cast<std::size_t>(n + 1), Dyad{0, 0});

    for (int x = 1; x <= n; ++x) {
        int len = bit_length(x);
        Dyad bestL{0, 0}, bestR{0, 0};
        bool hasL = false, hasR = false;

        for (int i = 0; i < len; ++i) {
            int bit = (x >> (len - 1 - i)) & 1;
            int y = static_cast<int>(delete_bit(x, len, i));
            const Dyad& vy = v[static_cast<std::size_t>(y)];

            if (bit == 1) {
                if (!hasL || cmp_dyad(vy, bestL) > 0) {
                    bestL = vy;
                    hasL = true;
                }
            } else {
                if (!hasR || cmp_dyad(vy, bestR) < 0) {
                    bestR = vy;
                    hasR = true;
                }
            }
        }

        v[static_cast<std::size_t>(x)] = simplest_between(bestL, hasL, bestR, hasR);
    }

    return v;
}

static i128 S_value(int n) {
    auto vals = compute_values(n);
    int max_exp = 0;
    for (int i = 1; i <= n; ++i) max_exp = std::max(max_exp, vals[static_cast<std::size_t>(i)].exp);

    i128 scaled_sum = 0;
    for (int i = 1; i <= n; ++i) {
        const auto& d = vals[static_cast<std::size_t>(i)];
        int sh = max_exp - d.exp;
        i128 term = static_cast<i128>(i) * static_cast<i128>(d.num);
        scaled_sum += (term << sh);
    }

    i128 den = static_cast<i128>(1) << max_exp;
    return (scaled_sum + den - 1) / den;
}

static std::string to_string_i128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    assert(S_value(2) == 2);
    assert(S_value(5) == 17);
    assert(S_value(10) == 64);

    std::cout << to_string_i128(S_value(100'000)) << '\n';
    return 0;
}
