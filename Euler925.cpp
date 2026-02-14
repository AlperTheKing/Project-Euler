#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kInv6 = 166666668ULL;

u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= kMod) a -= kMod;
    return a;
}

u64 sub_mod(u64 a, u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

u64 mul_mod(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

bool next_perm_digits(std::uint8_t* digits, int length) {
    int i = length - 2;
    while (i >= 0 && digits[i] >= digits[i + 1]) --i;
    if (i < 0) return false;
    int j = length - 1;
    while (digits[j] <= digits[i]) --j;
    std::swap(digits[i], digits[j]);
    std::reverse(digits + i + 1, digits + length);
    return true;
}

u64 parse_mod_from_digits(const std::uint8_t* digits, int length) {
    u64 v = 0;
    for (int i = 0; i < length; ++i) {
        v = (v * 10 + digits[i]) % kMod;
    }
    return v;
}

u64 b_diff_mod(u128 n) {
    if (n == 0) return 0;

    u128 x = n;
    std::uint8_t rev[64];
    int length = 0;
    while (x > 0) {
        rev[length++] = static_cast<std::uint8_t>(x % 10);
        x /= 10;
    }

    std::uint8_t digits[64];
    for (int i = 0; i < length; ++i) digits[i] = rev[length - 1 - i];

    const u64 n_mod = static_cast<u64>(n % kMod);
    if (!next_perm_digits(digits, length)) {
        return (n_mod == 0) ? 0 : (kMod - n_mod);
    }

    const u64 nxt_mod = parse_mod_from_digits(digits, length);
    return sub_mod(nxt_mod, n_mod);
}

u64 b_of_square_mod(u64 n) {
    const u128 sq = static_cast<u128>(n) * n;

    std::uint8_t rev[64];
    int length = 0;
    u128 x = sq;
    while (x > 0) {
        rev[length++] = static_cast<std::uint8_t>(x % 10);
        x /= 10;
    }

    if (length == 0) return 0;

    std::uint8_t digits[64];
    for (int i = 0; i < length; ++i) digits[i] = rev[length - 1 - i];
    if (!next_perm_digits(digits, length)) return 0;

    return parse_mod_from_digits(digits, length);
}

u64 brute_sum(u64 n) {
    u64 total = 0;
    for (u64 i = 1; i <= n; ++i) total = add_mod(total, b_of_square_mod(i));
    return total;
}

struct Context {
    int k;
    std::vector<u64> pow10_u64;
    std::vector<u128> pow10_u128;
    std::vector<u64> pow10_mod;
};

u64 monotone_diff(u64 root, int length, int trailing, const Context& ctx) {
    const int limit = ctx.k;
    const u128 p10t = ctx.pow10_u128[trailing];
    const u128 root_t = static_cast<u128>(root) * p10t;

    if (length + trailing >= limit) {
        return b_diff_mod(root_t * root_t);
    }

    const u128 p10t2 = p10t * p10t;
    const u64 p10a = ctx.pow10_u64[length];
    const u64 p10b = 10ULL * p10a;
    const u128 p10c = ctx.pow10_u128[limit - length - trailing - 1];
    const u128 p10d = static_cast<u128>(p10b) * p10t2;

    const u64 msb1 = static_cast<u64>((static_cast<u128>(root) * root) % p10a) / (p10a / 10ULL);

    u64 total = 0;
    for (int d = 0; d <= 9; ++d) {
        const u64 candidate = p10a * static_cast<u64>(d) + root;
        const u128 cand_sq = static_cast<u128>(candidate) * candidate;
        const u64 square = static_cast<u64>(cand_sq % p10b);
        const u128 square_t = static_cast<u128>(square) * p10t2;
        const u64 msb2 = square / p10a;

        if (msb2 < msb1) {
            const u64 diff_hi = b_diff_mod(p10d + square_t);
            if (cand_sq == square) {
                const u64 diff_lo = b_diff_mod(square_t);
                const u64 mul = static_cast<u64>((p10c - 1) % kMod);
                total = add_mod(total, diff_lo);
                total = add_mod(total, mul_mod(mul, diff_hi));
            } else {
                const u64 mul = static_cast<u64>(p10c % kMod);
                total = add_mod(total, mul_mod(mul, diff_hi));
            }
        } else {
            total = add_mod(total, monotone_diff(candidate, length + 1, trailing, ctx));
        }
    }

    return total;
}

u64 solve_925(int k) {
    Context ctx;
    ctx.k = k;
    ctx.pow10_u64.assign(k + 1, 1);
    ctx.pow10_u128.assign(2 * k + 2, 1);
    ctx.pow10_mod.assign(k + 1, 1);

    for (int i = 1; i <= k; ++i) {
        ctx.pow10_u64[i] = ctx.pow10_u64[i - 1] * 10ULL;
        ctx.pow10_mod[i] = mul_mod(ctx.pow10_mod[i - 1], 10ULL);
    }
    for (int i = 1; i <= 2 * k + 1; ++i) {
        ctx.pow10_u128[i] = ctx.pow10_u128[i - 1] * static_cast<u128>(10);
    }

    const u64 n_mod = ctx.pow10_mod[k];
    const u64 n1_mod = (n_mod + kMod - 1) % kMod;
    const u64 two_n1_mod = (2ULL * n_mod + kMod - 1) % kMod;

    u64 total = mul_mod(mul_mod(mul_mod(n_mod, n1_mod), two_n1_mod), kInv6);

    for (int d = 1; d <= 9; ++d) {
        for (int t = 0; t < k; ++t) {
            total = add_mod(total, monotone_diff(static_cast<u64>(d), 1, t, ctx));
        }
    }

    return total;
}

void validate() {
    assert(brute_sum(10) == 270);
    assert(brute_sum(100) == 335316);
    assert(solve_925(1) == 270);
    assert(solve_925(2) == 335316);
    assert(solve_925(6) == 265829902ULL);
}

}  // namespace

int main(int argc, char** argv) {
    int k = 16;
    if (argc > 1) k = std::atoi(argv[1]);
    validate();
    std::cout << solve_925(k) << '\n';
    return 0;
}
