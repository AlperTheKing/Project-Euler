#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kP = 1'000'000'007ULL;

static u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 res = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = static_cast<u64>((static_cast<u128>(res) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return res;
}

static bool contains_digit(u128 x, u64 base, u64 digit) {
    while (x > 0) {
        if (static_cast<u64>(x % base) == digit) {
            return true;
        }
        x /= base;
    }
    return digit == 0;
}

static u64 brute_M(u64 n, u64 d) {
    for (u64 m = 1;; ++m) {
        const u128 sq = static_cast<u128>(m) * m;
        if (contains_digit(sq, n, d)) {
            return m;
        }
    }
}

static u128 ceil_isqrt(u128 n) {
    if (n == 0) {
        return 0;
    }
    long double approx = std::sqrt(static_cast<long double>(n));
    u128 r = static_cast<u128>(approx);
    while (r * r < n) {
        ++r;
    }
    while (r > 0 && (r - 1) * (r - 1) >= n) {
        --r;
    }
    return r;
}

static u64 M_for_target_digit(u64 d, u64 exp_leg, u64 exp_sqrt, u128 p2) {
    const u64 a = kP - d;
    const u64 leg = mod_pow(a, exp_leg, kP);

    if (leg == 1) {
        const u64 r = mod_pow(a, exp_sqrt, kP);
        return std::min(r, kP - r);
    }

    u128 L = p2 - static_cast<u128>(d) * kP;
    const u128 width = kP;
    while (true) {
        const u128 m = ceil_isqrt(L);
        const u128 sq = m * m;
        if (sq <= L + (width - 1)) {
            return static_cast<u64>(m);
        }
        L += p2;
    }
}

static std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    assert(brute_M(10, 7) == 24);
    assert(brute_M(11, 10) == 19);

    const u64 exp_leg = (kP - 1) / 2;
    const u64 exp_sqrt = (kP + 1) / 4;
    const u128 p2 = static_cast<u128>(kP) * kP;

    u128 total = 0;
    for (u64 d = 1; d <= 100'000; ++d) {
        total += M_for_target_digit(d, exp_leg, exp_sqrt, p2);
    }

    std::cout << to_string_u128(total) << '\n';
    return 0;
}
