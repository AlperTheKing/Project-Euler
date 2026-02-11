#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1'000'000'007LL;

i64 mod_pow(i64 base, i64 exp) {
    i64 result = 1;
    while (exp > 0) {
        if (exp & 1LL) {
            result = static_cast<i64>((__int128)result * base % kMod);
        }
        base = static_cast<i64>((__int128)base * base % kMod);
        exp >>= 1;
    }
    return result;
}

u64 brute_sum(const int n_power) {
    const int limit = 1 << n_power;

    std::vector<std::uint8_t> a(static_cast<std::size_t>(limit + 1), 0);
    std::vector<std::uint8_t> b(static_cast<std::size_t>(limit + 1), 0);
    std::vector<std::uint8_t> c(static_cast<std::size_t>(limit + 1), 0);
    std::vector<std::uint8_t> d(static_cast<std::size_t>(limit + 1), 0);

    a[0] = 0;
    b[0] = 1;
    c[0] = 0;
    d[0] = 1;

    for (int r = 1; r <= limit; ++r) {
        std::uint8_t ar = 0;
        std::uint8_t br = 0;
        std::uint8_t cr = 1;
        std::uint8_t dr = 1;

        for (int x = 1; x <= r; ++x) {
            const int parity = std::popcount(static_cast<unsigned>(x)) & 1;
            if (parity == 0) {
                ar = static_cast<std::uint8_t>(ar | c[static_cast<std::size_t>(r - x)]);
                br = static_cast<std::uint8_t>(br | d[static_cast<std::size_t>(r - x)]);
                cr = static_cast<std::uint8_t>(cr & a[static_cast<std::size_t>(r - x)]);
                dr = static_cast<std::uint8_t>(dr & b[static_cast<std::size_t>(r - x)]);
            } else {
                ar = static_cast<std::uint8_t>(ar | d[static_cast<std::size_t>(r - x)]);
                br = static_cast<std::uint8_t>(br | c[static_cast<std::size_t>(r - x)]);
                cr = static_cast<std::uint8_t>(cr & b[static_cast<std::size_t>(r - x)]);
                dr = static_cast<std::uint8_t>(dr & a[static_cast<std::size_t>(r - x)]);
            }
        }

        a[static_cast<std::size_t>(r)] = ar;
        b[static_cast<std::size_t>(r)] = br;
        c[static_cast<std::size_t>(r)] = cr;
        d[static_cast<std::size_t>(r)] = dr;
    }

    u64 total = 0;
    for (int n = 1; n <= limit; ++n) {
        const int parity = std::popcount(static_cast<unsigned>(n)) & 1;
        const bool oscar_wins = parity == 0 ? (a[static_cast<std::size_t>(n)] != 0)
                                            : (b[static_cast<std::size_t>(n)] != 0);
        if (!oscar_wins) {
            total += static_cast<u64>(n);
        }
    }
    return total;
}

i64 solve_mod(const int n_power) {
    i64 ans = 0;

    i64 c_k = 1;     // 2^k
    i64 s_k = 0;     // sum(E_k)
    i64 pow4_k = 1;  // 4^k
    i64 odd_pow4 = 4;

    for (int m = 0; m < n_power; ++m) {
        if ((m & 1) == 0) {
            ans += static_cast<i64>(((__int128)c_k * pow4_k + s_k) % kMod);
            ans %= kMod;

            const i64 next_s = static_cast<i64>((2 * s_k + (__int128)(c_k + 2) * pow4_k) % kMod);
            s_k = next_s;
            c_k = static_cast<i64>((2 * c_k) % kMod);
            pow4_k = static_cast<i64>((4 * pow4_k) % kMod);
        } else {
            ans += odd_pow4 - 1;
            ans %= kMod;
            if (ans < 0) {
                ans += kMod;
            }
            odd_pow4 = static_cast<i64>((4 * odd_pow4) % kMod);
        }
    }

    if ((n_power & 1) == 0) {
        ans += mod_pow(2, n_power);
        ans %= kMod;
    }

    return ans;
}

}  // namespace

int main() {
    assert(brute_sum(4) == 46);
    assert(brute_sum(12) == 54'532);
    assert(solve_mod(12) == 54'532);

    std::cout << solve_mod(12'345'678) << '\n';
    return 0;
}
