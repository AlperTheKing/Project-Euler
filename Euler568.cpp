#include <cassert>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

using u64 = std::uint64_t;

struct dd {
    double hi;
    double lo;
};

// Error-free transforms for double-double arithmetic.
static inline dd two_sum(double a, double b) {
    double s = a + b;
    double bb = s - a;
    double err = (a - (s - bb)) + (b - bb);
    return {s, err};
}

static inline dd quick_two_sum(double a, double b) {
    // Assumes |a| >= |b|.
    double s = a + b;
    double err = b - (s - a);
    return {s, err};
}

static inline dd two_prod(double a, double b) {
    double p = a * b;
    // Dekker split for IEEE-754 double.
    static constexpr double split = 134217729.0;  // 2^27 + 1
    double ca = split * a;
    double a_hi = ca - (ca - a);
    double a_lo = a - a_hi;
    double cb = split * b;
    double b_hi = cb - (cb - b);
    double b_lo = b - b_hi;
    double err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    return {p, err};
}

static inline dd dd_add(dd x, dd y) {
    dd s = two_sum(x.hi, y.hi);
    double e = x.lo + y.lo;
    dd t = two_sum(s.lo, e);
    dd u = quick_two_sum(s.hi, t.hi);
    double lo = u.lo + t.lo;
    return quick_two_sum(u.hi, lo);
}

static inline dd dd_sub(dd x, dd y) { return dd_add(x, {-y.hi, -y.lo}); }

static inline dd dd_mul_d(dd x, double y) {
    dd p = two_prod(x.hi, y);
    double e = x.lo * y;
    dd s = two_sum(p.lo, e);
    dd u = quick_two_sum(p.hi, s.hi);
    double lo = u.lo + s.lo;
    return quick_two_sum(u.hi, lo);
}

static inline dd to_dd(long double x) {
    double hi = static_cast<double>(x);
    long double rem = x - static_cast<long double>(hi);
    double lo = static_cast<double>(rem);
    return {hi, lo};
}

static long double harmonic_asymptotic(u64 n) {
    // Euler–Maclaurin:
    // H_n = ln(n) + gamma + 1/(2n) - 1/(12n^2) + 1/(120n^4) - 1/(252n^6) + ...
    static constexpr long double gamma = 0.57721566490153286060651209008240243104215933593992L;
    long double nn = static_cast<long double>(n);
    long double inv = 1.0L / nn;
    long double inv2 = inv * inv;
    long double inv4 = inv2 * inv2;
    long double inv6 = inv4 * inv2;
    long double inv8 = inv4 * inv4;
    long double inv10 = inv8 * inv2;
    return logl(nn) + gamma + 0.5L * inv - (1.0L / 12.0L) * inv2 + (1.0L / 120.0L) * inv4 -
           (1.0L / 252.0L) * inv6 + (1.0L / 240.0L) * inv8 - (5.0L / 660.0L) * inv10;
}

static u64 leading7_D(u64 n) {
    // D(n) = J_B(n) - J_A(n) = H_n / 2^n.
    long double Hn = 0;
    if (n <= 2000000ULL) {
        for (u64 k = 1; k <= n; ++k) Hn += 1.0L / static_cast<long double>(k);
    } else {
        Hn = harmonic_asymptotic(n);
    }

    // log10(2) split into hi+lo so n*log10(2) keeps enough fractional precision.
    static constexpr dd log10_2 = {0.3010299956639812, -4.786261105275507e-18};

    dd log10_H = to_dd(log10l(Hn));
    dd nlog10_2 = dd_mul_d(log10_2, static_cast<double>(n));
    dd L = dd_sub(log10_H, nlog10_2);

    // frac(L) in [0,1).
    double frac = L.hi - std::floor(L.hi);
    frac += L.lo;
    frac -= std::floor(frac);

    long double x = powl(10.0L, static_cast<long double>(frac) + 6.0L);
    long double eps =
        8.0L * std::numeric_limits<long double>::epsilon() * fabsl(x);  // n=6 needs ~1e-8.
    u64 ans = static_cast<u64>(floorl(x + eps));
    if (ans >= 10000000ULL) ans = 9999999ULL;  // mantissa is in [1,10)
    return ans;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Statement example: D(6)=0.03828125 -> leading 7 digits are 3828125.
    assert(leading7_D(6) == 3828125ULL);

    std::cout << leading7_D(123456789ULL) << '\n';
    return 0;
}
