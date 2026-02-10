#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

// Project Euler 583: Heron Envelopes
//
// Let the envelope be the pentagon A-B-C-D-E where ABDE is a rectangle and BCD is an isosceles
// flap on top of the rectangle.
//
// Write the rectangle width as BD = 2a and height as AB = h (a,h > 0). The flap sides are
// BC = CD = s and its perpendicular height above BD is t.
//
// Integrality conditions:
// - Rectangle diagonals AD and BE are integral  <=>  (2a)^2 + h^2 is a square.
// - Diagonals AC and CE are integral.
//   With coordinates A(0,0), B(0,h), D(2a,h), C(a,h+t):
//     AC^2 = a^2 + (h+t)^2.
//   Hence AC integral  <=>  a^2 + (h+t)^2 is a square.
// - The flap itself gives the right triangle with legs a and t:
//     s^2 = a^2 + t^2.
//
// Also, a sensible envelope requires t < h.
// The perimeter is P = AB+BC+CD+DE+EA = 2h + 2s + 2a = 2(a+h+s).
//
// Therefore we need all integer triples (a,h,t) with t<h such that
//   a^2 + t^2 is a square,
//   a^2 + (h+t)^2 is a square,
//   (2a)^2 + h^2 is a square,
// and a+h+s <= floor(p/2), where s = sqrt(a^2+t^2).
//
// We enumerate Pythagorean triples using Euclid's formula, generate:
// - all (a,t,s) with s <= L (L=floor(p/2)),
// - all (a,h) with h <= L and (2a,h) being legs of a Pythagorean triple,
// then for each matching a test the remaining diagonal condition by a square check.

using u64 = std::uint64_t;

struct Tri {
    int a;
    int t;
    int s;
};

struct Rect {
    int a;
    int h;
};

static inline bool is_square_u64(u64 x) {
    u64 r = (u64)std::llround(std::sqrt((long double)x));
    while (r * r > x) --r;
    while ((r + 1) * (r + 1) <= x) ++r;
    return r * r == x;
}

static std::vector<Tri> generate_triangles(int L) {
    std::vector<Tri> tris;

    const int mmax = (int)std::sqrt((long double)L) + 2;
    for (int m = 2; m <= mmax; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0) continue;
            if (std::gcd(m, n) != 1) continue;

            const int a0 = m * m - n * n;
            const int b0 = 2 * m * n;
            const int c0 = m * m + n * n;
            if (c0 > L) continue;

            const int kmax = L / c0;
            for (int k = 1; k <= kmax; ++k) {
                const int a = k * a0;
                const int t = k * b0;
                const int s = k * c0;
                tris.push_back(Tri{a, t, s});
                tris.push_back(Tri{t, a, s});
            }
        }
    }

    std::sort(tris.begin(), tris.end(), [](const Tri& x, const Tri& y) {
        if (x.a != y.a) return x.a < y.a;
        return x.t < y.t;
    });
    return tris;
}

static std::vector<Rect> generate_rectangles(int L) {
    // Generate Pythagorean triples with legs <= 2L and keep those where one leg is 2a and the other is h<=L.
    const int max_leg = 2 * L;
    std::vector<Rect> rects;

    const int mmax = (int)std::sqrt((long double)max_leg) + 2;
    for (int m = 2; m <= mmax; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0) continue;
            if (std::gcd(m, n) != 1) continue;

            const int a0 = m * m - n * n;
            const int b0 = 2 * m * n;
            const int max0 = std::max(a0, b0);
            if (max0 > max_leg) continue;

            const int kmax = max_leg / max0;
            for (int k = 1; k <= kmax; ++k) {
                const int x = k * a0;
                const int y = k * b0;

                if ((x & 1) == 0) {
                    const int a = x / 2;
                    const int h = y;
                    if (a > 0 && a <= L && h > 0 && h <= L) rects.push_back(Rect{a, h});
                }
                if ((y & 1) == 0) {
                    const int a = y / 2;
                    const int h = x;
                    if (a > 0 && a <= L && h > 0 && h <= L) rects.push_back(Rect{a, h});
                }
            }
        }
    }

    std::sort(rects.begin(), rects.end(), [](const Rect& x, const Rect& y) {
        if (x.a != y.a) return x.a < y.a;
        return x.h < y.h;
    });
    return rects;
}

static u64 S(u64 p) {
    const int L = (int)(p / 2);
    const auto tris = generate_triangles(L);
    const auto rects = generate_rectangles(L);

    u64 sum = 0;

    size_t it = 0;
    size_t ir = 0;
    while (it < tris.size() && ir < rects.size()) {
        const int a_tri = tris[it].a;
        const int a_rec = rects[ir].a;
        if (a_tri < a_rec) {
            const int a = a_tri;
            while (it < tris.size() && tris[it].a == a) ++it;
            continue;
        }
        if (a_rec < a_tri) {
            const int a = a_rec;
            while (ir < rects.size() && rects[ir].a == a) ++ir;
            continue;
        }

        const int a = a_tri;
        const size_t it0 = it;
        while (it < tris.size() && tris[it].a == a) ++it;
        const size_t it1 = it;

        const size_t ir0 = ir;
        while (ir < rects.size() && rects[ir].a == a) ++ir;
        const size_t ir1 = ir;

        const u64 a2 = (u64)a * (u64)a;

        for (size_t r = ir0; r < ir1; ++r) {
            const int h = rects[r].h;
            if (h <= 0) continue;
            const int limit_s = L - a - h;
            if (limit_s <= 0) continue;

            for (size_t t = it0; t < it1; ++t) {
                const int tv = tris[t].t;
                if (tv >= h) break; // tris sorted by t
                const int s = tris[t].s;
                if (s > limit_s) continue;

                const int u = h + tv;
                const u64 val = a2 + (u64)u * (u64)u;
                if (!is_square_u64(val)) continue;

                sum += 2ULL * (u64)(a + h + s);
            }
        }
    }

    return sum;
}

int main() {
    // Validation point from the statement.
    const u64 s1e4 = S(10000);
    if (s1e4 != 884680ULL) {
        std::cerr << "Validation failed: S(1e4) got " << s1e4 << "\n";
        return 1;
    }

    std::cout << S(10000000ULL) << "\n";
    return 0;
}
