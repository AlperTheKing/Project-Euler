#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 398874989ULL;
constexpr u64 kGroup = kMod * kMod - 1ULL;

struct Elem {
    u64 a;
    u64 b;
};

u64 add_mod(u64 x, u64 y, u64 mod) {
    x += y;
    if (x >= mod) {
        x -= mod;
    }
    return x;
}

u64 mul_mod(u64 x, u64 y, u64 mod) {
    return static_cast<u64>((static_cast<u128>(x) * y) % mod);
}

Elem mul_elem(const Elem& x, const Elem& y) {
    const u64 aa = mul_mod(x.a, y.a, kMod);
    const u64 bb = mul_mod(x.b, y.b, kMod);
    const u64 real = add_mod(aa, (5ULL * bb) % kMod, kMod);
    const u64 imag = add_mod(mul_mod(x.a, y.b, kMod), mul_mod(x.b, y.a, kMod), kMod);
    return {real, imag};
}

std::array<Elem, 64> precompute_powers(const Elem& base) {
    std::array<Elem, 64> powers{};
    powers[0] = base;
    for (int i = 1; i < 64; ++i) {
        powers[i] = mul_elem(powers[i - 1], powers[i - 1]);
    }
    return powers;
}

Elem pow_from_powers(u64 exp, const std::array<Elem, 64>& powers) {
    Elem result{1, 0};
    while (exp != 0) {
        const int bit = __builtin_ctzll(exp);
        result = mul_elem(result, powers[bit]);
        exp &= (exp - 1);
    }
    return result;
}

Elem pow_elem(Elem base, u64 exp) {
    Elem result{1, 0};
    while (exp != 0) {
        if (exp & 1ULL) {
            result = mul_elem(result, base);
        }
        exp >>= 1U;
        if (exp != 0) {
            base = mul_elem(base, base);
        }
    }
    return result;
}

u64 pow_mod(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 value = base % mod;
    while (exp != 0) {
        if (exp & 1ULL) {
            result = mul_mod(result, value, mod);
        }
        exp >>= 1U;
        if (exp != 0) {
            value = mul_mod(value, value, mod);
        }
    }
    return result;
}

u64 pow5_mod(u64 x) {
    const u64 x2 = mul_mod(x, x, kMod);
    const u64 x4 = mul_mod(x2, x2, kMod);
    return mul_mod(x4, x, kMod);
}

u64 s_from_exp(u64 exp, const std::array<Elem, 64>& base_powers) {
    const Elem value = pow_from_powers(exp, base_powers);
    const u64 p = value.b;
    const u64 q = (value.a == 0 ? 0 : kMod - value.a);
    return add_mod(pow5_mod(p), pow5_mod(q), kMod);
}

u64 solve() {
    constexpr int kM = 1'618'034;
    const Elem base{(kMod + kMod - 2) % kMod, 1};  // -2 + sqrt(5)
    const auto base_powers = precompute_powers(base);

    u64 exp_prev2 = 5;  // 5^{F_1}
    u64 exp_prev1 = 5;  // 5^{F_2}

    u64 answer = s_from_exp(exp_prev1, base_powers);  // i=2
    for (int i = 3; i <= kM; ++i) {
        const u64 current = mul_mod(exp_prev1, exp_prev2, kGroup);
        answer = add_mod(answer, s_from_exp(current, base_powers), kMod);
        exp_prev2 = exp_prev1;
        exp_prev1 = current;
    }
    return answer;
}

void validate() {
    const Elem base{(kMod + kMod - 2) % kMod, 1};
    const auto base_powers = precompute_powers(base);

    assert(s_from_exp(1, base_powers) == 33);
    assert(s_from_exp(5, base_powers) == 257933744);
    assert(s_from_exp(25, base_powers) == 26500067);

    Elem stepped = base;
    u64 exp = 1;
    for (int n = 0; n <= 20; ++n) {
        const Elem by_exp = pow_from_powers(exp, base_powers);
        assert(by_exp.a == stepped.a && by_exp.b == stepped.b);
        stepped = pow_elem(stepped, 5);
        exp = mul_mod(exp, 5, kGroup);
    }

    u64 fib_prev2 = 1;
    u64 fib_prev1 = 1;
    u64 fexp_prev2 = 5;
    u64 fexp_prev1 = 5;
    for (int i = 3; i <= 25; ++i) {
        const u64 fib_cur = fib_prev1 + fib_prev2;
        const u64 fexp_cur = mul_mod(fexp_prev1, fexp_prev2, kGroup);
        assert(fexp_cur == pow_mod(5, fib_cur, kGroup));
        fib_prev2 = fib_prev1;
        fib_prev1 = fib_cur;
        fexp_prev2 = fexp_prev1;
        fexp_prev1 = fexp_cur;
    }
}

}  // namespace

int main() {
    validate();
    std::cout << solve() << '\n';
    return 0;
}
