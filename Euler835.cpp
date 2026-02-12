#include <cassert>
#include <cstdint>
#include <iostream>

#include <boost/multiprecision/cpp_dec_float.hpp>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using Real = boost::multiprecision::cpp_dec_float_100;

static constexpr u64 kMod = 1'234'567'891ULL;

static u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) r = static_cast<u64>((static_cast<u128>(r) * a) % mod);
        a = static_cast<u64>((static_cast<u128>(a) * a) % mod);
        e >>= 1ULL;
    }
    return r;
}

static u64 sum_family_I_mod(u64 half_exp) {
    // N = 10^(2*half_exp), s = 10^half_exp, q = s/2
    const u64 s_mod = mod_pow(10, half_exp, kMod);
    const u64 inv2 = (kMod + 1) / 2;
    const u64 inv3 = mod_pow(3, kMod - 2, kMod);
    const u64 q = mod_mul(s_mod, inv2);

    // Sum over odd t = 1..s-1: t(t+1) = q(q+1)(4q-1)/3
    u64 part = q;
    part = mod_mul(part, (q + 1) % kMod);
    u64 four_q_minus_one = (mod_mul(4 % kMod, q) + kMod - 1) % kMod;
    part = mod_mul(part, four_q_minus_one);
    part = mod_mul(part, inv3);

    // remove t=1 term (degenerate triangle)
    part = (part + kMod - 2) % kMod;
    return part;
}

struct Mat3 {
    u64 a[3][3]{};
};

static Mat3 mat_mul(const Mat3& A, const Mat3& B) {
    Mat3 C{};
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            if (A.a[i][k] == 0) continue;
            for (int j = 0; j < 3; ++j) {
                if (B.a[k][j] == 0) continue;
                C.a[i][j] = (C.a[i][j] + static_cast<u64>((static_cast<u128>(A.a[i][k]) * B.a[k][j]) % kMod)) % kMod;
            }
        }
    }
    return C;
}

static Mat3 mat_pow(Mat3 base, u64 exp) {
    Mat3 res{};
    for (int i = 0; i < 3; ++i) res.a[i][i] = 1;
    while (exp > 0) {
        if (exp & 1ULL) res = mat_mul(res, base);
        base = mat_mul(base, base);
        exp >>= 1ULL;
    }
    return res;
}

static u64 sum_family_II_mod(u64 K) {
    if (K == 0) return 0;
    if (K == 1) return 12;

    // p_{k+1} = 6 p_k - p_{k-1}
    // s_{k+1} = s_k + p_{k+1}
    Mat3 A{};
    A.a[0][0] = 6;
    A.a[0][1] = kMod - 1;
    A.a[0][2] = 0;

    A.a[1][0] = 1;
    A.a[1][1] = 0;
    A.a[1][2] = 0;

    A.a[2][0] = 6;
    A.a[2][1] = kMod - 1;
    A.a[2][2] = 1;

    Mat3 P = mat_pow(A, K - 2);

    // state at k=2: [p2, p1, s2] = [70,12,82]
    const u64 v[3] = {70, 12, 82};
    u64 out[3] = {0, 0, 0};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out[i] = (out[i] + static_cast<u64>((static_cast<u128>(P.a[i][j]) * v[j]) % kMod)) % kMod;
        }
    }

    return out[2];
}

static u64 K_max_for_power10_N(u64 exp10) {
    // p_k = alpha * lambda^k + beta * mu^k, with lambda = 3 + 2*sqrt(2)
    const Real sqrt2 = sqrt(Real(2));
    const Real lambda = Real(3) + Real(2) * sqrt2;
    const Real mu = Real(3) - Real(2) * sqrt2;
    const Real p1 = 12;
    const Real p2 = 70;

    const Real alpha = (p2 - p1 * mu) / (lambda * (lambda - mu));

    const Real log10_lambda = log(lambda) / log(Real(10));
    const Real log10_alpha = log(alpha) / log(Real(10));

    const Real E = Real(exp10);
    Real kest = (E - log10_alpha) / log10_lambda;
    u64 K = kest.convert_to<u64>();

    auto log10_pk = [&](u64 k) -> Real {
        return Real(k) * log10_lambda + log10_alpha;
    };

    while (log10_pk(K + 1) <= E) ++K;
    while (K > 0 && log10_pk(K) > E) --K;

    return K;
}

static u64 S_bruteforce(u64 N) {
    u64 sum = 0;

    // Family I: (t, (t^2-1)/2, (t^2+1)/2), t odd >=3, perimeter t(t+1)
    for (u64 t = 3;; t += 2) {
        u128 p = static_cast<u128>(t) * (t + 1);
        if (p > N) break;
        sum += static_cast<u64>(p);
    }

    // Family II: consecutive legs sequence perimeters
    u64 u = 7, c = 5;  // first valid gives perimeter 12
    while (true) {
        u64 p = u + c;
        if (p > N) break;
        sum += p;
        u64 nu = 3 * u + 4 * c;
        u64 nc = 2 * u + 3 * c;
        u = nu;
        c = nc;
    }

    // overlap at 3-4-5 only
    sum -= 12;
    return sum;
}

int main() {
    assert(S_bruteforce(100) == 258);
    assert(S_bruteforce(10000) == 172004);

    const u64 exp10 = 10'000'000'000ULL;
    const u64 K = K_max_for_power10_N(exp10);

    const u64 sumI = sum_family_I_mod(exp10 / 2);
    const u64 sumII = sum_family_II_mod(K);

    u64 ans = (sumI + sumII) % kMod;
    ans = (ans + kMod - 12) % kMod;

    std::cout << ans << '\n';
    return 0;
}
