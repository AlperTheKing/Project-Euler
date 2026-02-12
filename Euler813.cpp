#include <cassert>
#include <cstdint>
#include <iostream>

#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::bit_test;
using boost::multiprecision::cpp_int;
using boost::multiprecision::lsb;
using boost::multiprecision::msb;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kMod = 1'000'000'007ULL;

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

static cpp_int clmul(cpp_int a, const cpp_int& b) {
    cpp_int res = 0;
    std::size_t shift = 0;
    while (a != 0) {
        const unsigned tz = static_cast<unsigned>(lsb(a));
        shift += tz;
        a >>= tz;
        res ^= (b << shift);
        a >>= 1;
        ++shift;
    }
    return res;
}

static cpp_int clpow(cpp_int base, u64 exp) {
    cpp_int res = 1;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = clmul(res, base);
        }
        exp >>= 1ULL;
        if (exp > 0) {
            base = clmul(base, base);
        }
    }
    return res;
}

static u64 eval_poly_mod(const cpp_int& poly_bits, u64 x_mod) {
    const std::size_t deg = static_cast<std::size_t>(msb(poly_bits));
    u64 ans = 0;
    for (std::size_t i = deg + 1; i-- > 0;) {
        ans = static_cast<u64>((static_cast<u128>(ans) * x_mod) % kMod);
        if (bit_test(poly_bits, i)) {
            ++ans;
            if (ans == kMod) {
                ans = 0;
            }
        }
    }
    return ans;
}

int main() {
    const cpp_int eleven = 11;
    assert(clmul(eleven, eleven) == 69);
    assert(clpow(eleven, 2) == 69);

    const u64 three_pow_8 = 6561ULL;
    const u64 two_pow_52 = (1ULL << 52);

    const cpp_int f = (cpp_int(1) << 3) ^ (cpp_int(1) << 1) ^ cpp_int(1);
    const cpp_int g = clpow(f, three_pow_8);

    const u64 x = mod_pow(2, two_pow_52);
    std::cout << eval_poly_mod(g, x) << '\n';
    return 0;
}
