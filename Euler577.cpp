#include <cstdint>
#include <iostream>
#include <string>

// Project Euler 577: Counting Hexagons
//
// Use axial coordinates (x,y) on the triangular lattice:
//   0 <= x <= n, 0 <= y <= n-x
// so there are (n+1)(n+2)/2 lattice points.
//
// A regular hexagon with vertices on this lattice is determined by a nonzero side vector
// u = (a,b) in the lattice (Eisenstein integers). Let v be u rotated by +60 degrees:
//   v = rot60(u) = (-b, a+b).
// Then for some base vertex P0 the hexagon vertices are:
//   P0,
//   P0 + u,
//   P0 + u + v,
//   P0 + 2v,
//   P0 + 2v - u,
//   P0 + v - u.
//
// To count each hexagon once, restrict to canonical u with a>0 and b>=0.
// For such u, the number of valid translations P0 that keep all vertices inside the big triangle
// depends only on s = a+b, and is:
//   placements(n,s) = (n-3s+1)(n-3s+2)/2  for 3s <= n, else 0.
// There are exactly s choices of (a,b) with a>0, b>=0 and a+b=s.
// Hence:
//   H(n) = sum_{s=1..floor(n/3)} s * (n-3s+1)(n-3s+2)/2.
//
// We need sum_{n=3..N} H(n). Swap sums with m = n-3s+1:
//   sum_{n=3..N} H(n) = sum_{s=1..floor(N/3)} s * sum_{m=1..(N-3s+1)} m(m+1)/2
//                    = sum_{s} s * M(M+1)(M+2)/6, where M = N-3s+1.

using i128 = __int128_t;

static std::string to_string_i128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x > 0) {
        int digit = int(x % 10);
        s.push_back(char('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static i128 H(int n) {
    i128 tot = 0;
    for (i128 s = 1; 3 * s <= n; ++s) {
        i128 m = (i128)n - 3 * s + 1;
        tot += s * m * (m + 1) / 2;
    }
    return tot;
}

static i128 sum_H_up_to(int N) {
    i128 ans = 0;
    for (i128 s = 1; 3 * s <= N; ++s) {
        i128 M = (i128)N - 3 * s + 1;
        ans += s * M * (M + 1) * (M + 2) / 6;
    }
    return ans;
}

int main() {
    // Validation from the statement.
    if (H(3) != 1 || H(6) != 12 || H(20) != 966) {
        std::cerr << "Validation failed\n";
        return 1;
    }

    static constexpr int N = 12345;
    std::cout << to_string_i128(sum_H_up_to(N)) << "\n";
    return 0;
}

