#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

// Project Euler 582: Integer-sided 120-degree triangles with b-a <= 100
//
// With a <= b <= c and the 120-degree angle opposite c, the law of cosines gives:
//   c^2 = a^2 + b^2 + ab.
//
// Parameterization via Eisenstein integers:
// Let m > n > 0 and define (a0,b0,c0):
//   a0 = m^2 - n^2
//   b0 = 2mn + n^2
//   c0 = m^2 + mn + n^2
// Then c0^2 = a0^2 + b0^2 + a0*b0.
// Every integer solution with gcd(a,b)=1 is a unit multiple of such a triple, and scaling by k
// gives all solutions:
//   (a,b,c) = k*(a0,b0,c0).
//
// The difference satisfies:
//   b0 - a0 = 2n^2 + 2mn - m^2.
// Setting X = m - n > 0, this becomes:
//   b0 - a0 = -(X^2 - 3n^2) = t,
// i.e. solutions correspond to Pell-type equations:
//   X^2 - 3n^2 = -t,  with 1 <= t <= 100.
//
// For fixed t, all positive solutions (X,n) are generated from finitely many fundamental ones
// by multiplication with the fundamental unit 2+sqrt(3):
//   (X',n') = (2X + 3n, X + 2n).
//
// For each solution we compute c0 and count admissible scaling factors k:
//   1 <= k <= floor(100/t) and k*c0 <= N  =>  k <= min(floor(100/t), floor(N/c0)).
//
// We sum this over all t and all solutions.

using boost::multiprecision::cpp_int;
using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Sol {
    u64 X = 0;
    u64 n = 0;
    bool operator<(const Sol& o) const { return (X < o.X) || (X == o.X && n < o.n); }
};

static Sol reduce_solution(Sol s) {
    // Multiply by (2 - sqrt(3)) repeatedly while staying positive to reach a minimal orbit rep.
    while (true) {
        const i64 X = (i64)s.X;
        const i64 n = (i64)s.n;
        const i64 Xp = 2 * X - 3 * n;
        const i64 np = -X + 2 * n;
        if (Xp <= 0 || np <= 0) break;
        s.X = (u64)Xp;
        s.n = (u64)np;
    }
    return s;
}

static std::vector<Sol> fundamental_solutions(int S) {
    // Brute small n to find solutions, then reduce to orbit representatives.
    std::set<Sol> reps;
    for (int n = 1; n <= 500; ++n) {
        const long long v = 3LL * n * n + S;
        if (v <= 0) continue;
        const long long x = (long long)std::llround(std::floor(std::sqrt((long double)v)));
        for (long long X = std::max(1LL, x - 2); X <= x + 2; ++X) {
            if (X * X != v) continue;
            Sol s{(u64)X, (u64)n};
            reps.insert(reduce_solution(s));
        }
    }
    return std::vector<Sol>(reps.begin(), reps.end());
}

static cpp_int c0_from_Xn(const cpp_int& X, const cpp_int& n) {
    // c0 = 3n^2 + 3Xn + X^2
    return 3 * n * n + 3 * X * n + X * X;
}

static u64 T(const cpp_int& N) {
    struct Tri {
        cpp_int a, b, c;
        bool operator<(const Tri& o) const {
            if (a != o.a) return a < o.a;
            if (b != o.b) return b < o.b;
            return c < o.c;
        }
    };

    std::set<Tri> uniq;

    // Cache fundamental solutions per signed RHS S (X^2 - 3n^2 = S).
    std::map<int, std::vector<Sol>> cache;

    for (int S = -100; S <= 100; ++S) {
        if (S == 0) continue;
        const int diff = std::abs(S);
        const int kmax = 100 / diff;
        if (kmax == 0) continue;

        auto& reps = cache[S];
        if (reps.empty()) reps = fundamental_solutions(S);

        for (const auto& rep : reps) {
            cpp_int X = rep.X;
            cpp_int n = rep.n;

            while (true) {
                const cpp_int m = n + X;

                const cpp_int A = m * m - n * n;             // m^2 - n^2
                const cpp_int B = 2 * m * n + n * n;         // 2mn + n^2
                const cpp_int C = m * m + m * n + n * n;     // m^2 + mn + n^2
                if (C > N) break;

                const cpp_int a0 = std::min(A, B);
                const cpp_int b0 = std::max(A, B);

                for (int k = 1; k <= kmax; ++k) {
                    const cpp_int kk = k;
                    const cpp_int c = kk * C;
                    if (c > N) break;
                    uniq.insert(Tri{kk * a0, kk * b0, c});
                }

                // Next (X,n): multiply by 2+sqrt(3).
                const cpp_int Xnxt = 2 * X + 3 * n;
                const cpp_int nnxt = X + 2 * n;
                X = Xnxt;
                n = nnxt;
            }
        }
    }

    return (u64)uniq.size();
}

static cpp_int pow10(int e) {
    cpp_int x = 1;
    for (int i = 0; i < e; ++i) x *= 10;
    return x;
}

int main() {
    // Validation points from the statement.
    const u64 t1000 = T(cpp_int(1000));
    if (t1000 != 235ULL) {
        std::cerr << "Validation failed: T(1000) got " << t1000 << "\n";
        return 1;
    }
    const u64 t1e8 = T(cpp_int(100000000));
    if (t1e8 != 1245ULL) {
        std::cerr << "Validation failed: T(1e8) got " << t1e8 << "\n";
        return 1;
    }

    const cpp_int N = pow10(100);
    std::cout << T(N) << "\n";
    return 0;
}
