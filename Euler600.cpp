#include <cstdint>
#include <algorithm>
#include <iostream>

// Project Euler 600: Integer Sided Equiangular Hexagons
//
// An equiangular hexagon has all interior angles 120 degrees, so successive edge directions
// differ by 60 degrees. With side lengths (a,b,c,d,e,f) around the hexagon, the closure
// condition becomes:
//   a - d = c - f,
//   b - e = f - c.
//
// We count ordered side-length 6-tuples of positive integers satisfying closure and
// perimeter <= n, then quotient by the dihedral symmetries of a hexagon (rotations and
// reflections) using Burnside's lemma.
//
// Burnside requires the number of tuples fixed by each symmetry. The rotation cases reduce
// to simple compositions. Reflections split into two types and reduce to tractable integer
// sums. The identity count can be computed by summing over the difference parameter p:
// writing d=A, f=C, e=E and letting
//   a = A + p,  b = E - p,  c = C + p,
// yields perimeter 2(A+C+E)+p and positivity bounds that depend only on p.

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i64 = std::int64_t;

static void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }
    char buf[64];
    int n = 0;
    while (x > 0) {
        buf[n++] = (char)('0' + (u64)(x % 10));
        x /= 10;
    }
    while (n--) std::cout << buf[n];
}

static inline u128 choose2(u64 t) {
    if (t < 2) return 0;
    return (u128)t * (u128)(t - 1) / 2;
}

static inline u128 choose3(u64 t) {
    if (t < 3) return 0;
    return (u128)t * (u128)(t - 1) * (u128)(t - 2) / 6;
}

static inline u128 tet(u64 r) {
    // Number of nonnegative integer triples with sum <= r is C(r+3,3).
    return (u128)(r + 1) * (u128)(r + 2) * (u128)(r + 3) / 6;
}

static u128 fix_identity(u64 n) {
    // Sum over p (called x in many derivations).
    // For fixed p, let A=d, C=f, E=e. Then:
    //   a=A+p >=1  => A >= 1-p
    //   c=C+p >=1  => C >= 1-p
    //   b=E-p >=1  => E >= p+1
    // Perimeter: 2(A+C+E)+p <= n  => A+C+E <= floor((n-p)/2).

    u128 total = 0;
    for (i64 p = -(i64)n; p <= (i64)n; ++p) {
        const i64 T = ((i64)n - p) / 2;
        const i64 A0 = std::max<i64>(1, 1 - p);
        const i64 E0 = std::max<i64>(1, p + 1);
        const i64 base = 2 * A0 + E0;
        if (T < base) continue;
        const i64 R = T - base;
        total += tet((u64)R);
    }
    return total;
}

static u128 fix_reflection_even(u64 n) {
    // Representative reflection constraints: (a,b,c,d,e,f) = (a,b,c,d,c,b).
    // Closure gives d = a - c + b.
    // Perimeter: 2a + 3b + c <= n, and d>=1 implies c <= a+b-1.

    u128 total = 0;
    if (n < 6) return 0;

    const u64 bmax = (n - 3) / 3; // need 3b+3 <= n for any solution
    for (u64 b = 1; b <= bmax; ++b) {
        const i64 M = (i64)n - 3 * (i64)b; // M = n - 3b

        // Need M - 2a >= 1 for c>=1.
        const i64 Amax = (M - 1) / 2;
        if (Amax <= 0) continue;

        // Split where a+b-1 <= M-2a  <=>  3a <= M-b+1.
        const i64 Asplit = ((i64)n - 4 * (i64)b + 1) / 3;
        i64 m = std::min<i64>(Amax, std::max<i64>(0, Asplit));

        // Region 1: a=1..m, c_max = a+b-1.
        if (m >= 1) {
            const u128 sum_a = (u128)m * (u128)(m + 1) / 2;
            total += sum_a + (u128)(b - 1) * (u128)m;
        } else {
            m = 0;
        }

        // Region 2: a=m+1..Amax, c_max = M-2a.
        if (Amax > m) {
            const i64 cnt = Amax - m;
            const u128 sum_a = (u128)Amax * (u128)(Amax + 1) / 2 - (u128)m * (u128)(m + 1) / 2;
            total += (u128)M * (u128)cnt - (u128)2 * sum_a;
        }
    }
    return total;
}

static u128 fix_reflection_odd(u64 n) {
    // Representative odd reflection fixes: (a,b,c,d,e,f)=(a,a,c,a,a,c).
    // Condition: 4a+2c <= n  <=> 2a + c <= floor(n/2).

    const u64 T = n / 2;
    if (T < 3) return 0;
    const u64 amax = (T - 1) / 2;
    // Sum_{a=1..amax} (T - 2a)
    return (u128)amax * (u128)(T - (amax + 1));
}

static u128 H(u64 n) {
    const u128 id = fix_identity(n);
    const u128 r1 = (u64)(n / 6);
    const u64 T3 = n / 3;
    const u128 r2 = choose2(T3);
    const u64 T2 = n / 2;
    const u128 r3 = choose3(T2);
    const u128 fe = fix_reflection_even(n);
    const u128 fo = fix_reflection_odd(n);

    const u128 num = id + 2 * r1 + 2 * r2 + r3 + 3 * fe + 3 * fo;
    return num / 12;
}

int main() {
    // Validation points from the statement.
    if (H(6) != 1) {
        std::cerr << "Validation failed: H(6)\n";
        return 1;
    }
    if (H(12) != 10) {
        std::cerr << "Validation failed: H(12)\n";
        return 1;
    }
    if (H(100) != 31248) {
        std::cerr << "Validation failed: H(100)\n";
        return 1;
    }

    print_u128(H(55106));
    std::cout << "\n";
    return 0;
}
