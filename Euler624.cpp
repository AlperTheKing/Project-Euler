#include <cassert>
#include <cstdint>
#include <iostream>

using u64 = unsigned long long;
using u128 = __uint128_t;

struct Mat2 {
    u64 a00, a01, a10, a11;
};

static u64 mod_mul(u64 a, u64 b, u64 mod) { return (u128)a * b % mod; }

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    while (e) {
        if (e & 1) r = mod_mul(r, a, mod);
        a = mod_mul(a, a, mod);
        e >>= 1;
    }
    return r;
}

static Mat2 mat_mul(const Mat2 &x, const Mat2 &y, u64 mod) {
    Mat2 r;
    r.a00 = (mod_mul(x.a00, y.a00, mod) + mod_mul(x.a01, y.a10, mod)) % mod;
    r.a01 = (mod_mul(x.a00, y.a01, mod) + mod_mul(x.a01, y.a11, mod)) % mod;
    r.a10 = (mod_mul(x.a10, y.a00, mod) + mod_mul(x.a11, y.a10, mod)) % mod;
    r.a11 = (mod_mul(x.a10, y.a01, mod) + mod_mul(x.a11, y.a11, mod)) % mod;
    return r;
}

static Mat2 mat_pow(Mat2 base, u64 e, u64 mod) {
    Mat2 r{1, 0, 0, 1};
    while (e) {
        if (e & 1) r = mat_mul(r, base, mod);
        base = mat_mul(base, base, mod);
        e >>= 1;
    }
    return r;
}

static u64 P_mod(u64 n, u64 mod) {
    assert(mod > 2);
    const u64 inv2 = (mod + 1) / 2;
    const Mat2 A{inv2, inv2, inv2, 0};

    const Mat2 An = mat_pow(A, n, mod);
    const Mat2 An1 = mat_pow(A, n - 1, mod);

    const u64 w0 = An1.a00;
    const u64 w1 = An1.a10;

    const u64 m00 = (1 + mod - An.a00) % mod;
    const u64 m01 = (0 + mod - An.a01) % mod;
    const u64 m10 = (0 + mod - An.a10) % mod;
    const u64 m11 = (1 + mod - An.a11) % mod;

    const u64 det = (mod_mul(m00, m11, mod) + mod - mod_mul(m01, m10, mod)) % mod;
    const u64 inv_det = mod_pow(det, mod - 2, mod);

    const u64 s1 = (mod_mul((mod - m10) % mod, w0, mod) + mod_mul(m00, w1, mod)) % mod;
    return mod_mul(inv2, mod_mul(s1, inv_det, mod), mod);
}

static u64 Q_P(u64 n, u64 mod) {
    u64 r = P_mod(n, mod);
    return (r == 0) ? mod : r;
}

int main() {
    assert(Q_P(2, 109) == 66);
    assert(Q_P(3, 109) == 46);
    std::cout << Q_P(1'000'000'000'000'000'000ULL, 1'000'000'009ULL) << "\n";
    return 0;
}

