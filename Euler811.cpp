#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kMod = 1'000'062'031ULL;

static u64 mod_pow(u64 base, u64 exp) {
    u64 res = 1;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = static_cast<u64>((static_cast<u128>(res) * base) % kMod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % kMod);
        exp >>= 1ULL;
    }
    return res;
}

static int bit_length(u64 x) {
    if (x == 0) {
        return 0;
    }
    return 64 - __builtin_clzll(x);
}

static std::vector<u64> binom_row(int r) {
    std::vector<u64> c(r + 1, 0);
    c[0] = 1;
    for (int k = 1; k <= r; ++k) {
        c[k] = static_cast<u64>((static_cast<u128>(c[k - 1]) * static_cast<u64>(r - k + 1)) /
                                static_cast<u64>(k));
    }
    return c;
}

static std::vector<u64> bit_positions_desc_u64(u64 n) {
    std::vector<u64> pos;
    for (int b = 63; b >= 0; --b) {
        if ((n >> b) & 1ULL) {
            pos.push_back(static_cast<u64>(b));
        }
    }
    return pos;
}

static u64 a_from_positions_desc(const std::vector<u64>& pos_desc) {
    const std::size_t ones = pos_desc.size();
    std::vector<u64> v(ones + 1, 0);
    v[0] = 1;
    for (std::size_t i = 1; i <= ones; ++i) {
        v[i] = (5ULL * v[i - 1] + 3ULL) % kMod;
    }

    u64 ans = 1;
    for (std::size_t i = 0; i < ones; ++i) {
        const u64 gap = (i + 1 < ones) ? (pos_desc[i] - pos_desc[i + 1] - 1ULL) : pos_desc[i];
        if (gap > 0) {
            ans = static_cast<u64>((static_cast<u128>(ans) * mod_pow(v[i + 1], gap)) % kMod);
        }
    }
    return ans;
}

static u64 a_bruteforce_mod(u64 n) {
    std::vector<u64> a(n + 1, 0);
    a[0] = 1;
    for (u64 i = 1; i <= n; ++i) {
        if (i & 1ULL) {
            a[i] = a[i >> 1ULL];
        } else {
            const u64 m = i >> 1ULL;
            const u64 low = m & (~m + 1ULL);
            a[i] = (3ULL * a[m] + 5ULL * a[i - low]) % kMod;
        }
    }
    return a[n];
}

static u64 h_no_overlap(u64 t, int r) {
    const std::vector<u64> coeff = binom_row(r);
    int max_coeff_bits = 0;
    for (u64 c : coeff) {
        max_coeff_bits = std::max(max_coeff_bits, bit_length(c));
    }
    assert(t > static_cast<u64>(max_coeff_bits));

    std::vector<u64> pos;
    for (int k = 0; k <= r; ++k) {
        u64 x = coeff[k];
        while (x > 0) {
            const int b = __builtin_ctzll(x);
            pos.push_back(static_cast<u64>(k) * t + static_cast<u64>(b));
            x &= (x - 1ULL);
        }
    }

    std::sort(pos.begin(), pos.end(), std::greater<u64>());
    return a_from_positions_desc(pos);
}

int main() {
    assert((24ULL & (~24ULL + 1ULL)) == 8ULL);

    for (u64 n = 0; n <= 5000; ++n) {
        const auto pos = bit_positions_desc_u64(n);
        assert(a_from_positions_desc(pos) == a_bruteforce_mod(n));
    }

    assert(h_no_overlap(3, 2) == 636056ULL);

    constexpr u64 t = 100'000'000'000'031ULL;
    constexpr int r = 62;
    std::cout << h_no_overlap(t, r) << '\n';
    return 0;
}
