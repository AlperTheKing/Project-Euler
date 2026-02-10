#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

// Project Euler 595: Incremental Random Sort
//
// After each check, any adjacent cards that form an increasing consecutive run (i, i+1, i+2, ...)
// are glued into a single bundle. A shuffle then produces a uniformly random permutation of the
// current bundles (bundle sizes do not matter for the randomization).
//
// Each bundle is always an interval of consecutive labels, in increasing order.
// If there are m bundles, label them 1..m in increasing-value order of those intervals.
// A shuffle is then a uniformly random permutation of 1..m.
// The next gluing step merges exactly the adjacent pairs (i, i+1) that appear consecutively as
// i immediately followed by i+1 in that permutation.
// Hence the state is fully determined by the bundle count m.
//
// Let P(m,b) be the probability that a random permutation of 1..m, after gluing consecutive
// increasing adjacencies, results in exactly b bundles.
// A permutation yields b bundles iff it can be seen as a permutation of b maximal consecutive
// blocks (a composition of m into b parts) with the constraint that block i is not immediately
// followed by block i+1 for any i.
// This factors as:
//   count(m,b) = C(m-1, b-1) * f(b),
// where f(b) is the number of permutations of 1..b with no adjacency i immediately followed by i+1.
// By inclusion-exclusion on the (b-1) forbidden adjacencies:
//   f(b) = sum_{r=0}^{b-1} (-1)^r * C(b-1, r) * (b-r)!.
//
// Let E_m be the expected number of shuffles needed to finish when currently at m bundles
// (after a gluing step). Then E_1 = 0 and for m>1:
//   E_m = 1 + sum_b P(m,b) * E_b,
// which rearranges to:
//   E_m = (1 + sum_{b<m} P(m,b) E_b) / (1 - P(m,m)).
//
// The initial deck is a random permutation of n singletons, so after the initial gluing the
// bundle count has distribution P(n,·). Therefore S(n) = sum_b P(n,b) E_b.

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_dec_float_100;

static cpp_int factorial(int n) {
    cpp_int r = 1;
    for (int i = 2; i <= n; ++i) r *= i;
    return r;
}

int main() {
    constexpr int N = 52;

    std::vector<cpp_int> fact(N + 1);
    fact[0] = 1;
    for (int i = 1; i <= N; ++i) fact[i] = fact[i - 1] * i;

    auto comb = [&](int n, int k) -> cpp_int {
        if (k < 0 || k > n) return 0;
        return fact[n] / (fact[k] * fact[n - k]);
    };

    // f(b): permutations of 1..b with no adjacent pattern i,i+1.
    std::vector<cpp_int> f(N + 1);
    f[0] = 1;
    for (int b = 1; b <= N; ++b) {
        cpp_int fb = 0;
        for (int r = 0; r <= b - 1; ++r) {
            const cpp_int term = comb(b - 1, r) * fact[b - r];
            if (r & 1)
                fb -= term;
            else
                fb += term;
        }
        f[b] = fb;
    }

    // Transition probabilities P(m,b).
    std::vector<std::vector<cpp_dec_float_100>> P(N + 1);
    for (int m = 1; m <= N; ++m) {
        P[m].assign(m + 1, cpp_dec_float_100(0));
        const cpp_dec_float_100 denom = cpp_dec_float_100(fact[m]);
        for (int b = 1; b <= m; ++b) {
            const cpp_int cnt = comb(m - 1, b - 1) * f[b];
            P[m][b] = cpp_dec_float_100(cnt) / denom;
        }
    }

    // Expected shuffles from each state.
    std::vector<cpp_dec_float_100> E(N + 1);
    E[1] = 0;
    for (int m = 2; m <= N; ++m) {
        cpp_dec_float_100 num = 1;
        for (int b = 1; b < m; ++b) num += P[m][b] * E[b];
        E[m] = num / (cpp_dec_float_100(1) - P[m][m]);
    }

    auto S = [&](int n) -> cpp_dec_float_100 {
        cpp_dec_float_100 s = 0;
        for (int b = 1; b <= n; ++b) s += P[n][b] * E[b];
        return s;
    };

    // Validation points from the statement.
    {
        const auto s1 = S(1);
        if (s1 != 0) {
            std::cerr << "Validation failed: S(1)\n";
            return 1;
        }
        const auto s2 = S(2);
        if (s2 != 1) {
            std::cerr << "Validation failed: S(2)\n";
            return 1;
        }
        const auto s5 = S(5);
        const cpp_dec_float_100 target = cpp_dec_float_100(4213) / cpp_dec_float_100(871);
        if (abs(s5 - target) > cpp_dec_float_100("1e-60")) {
            std::cerr << "Validation failed: S(5)\n";
            return 1;
        }
    }

    const auto ans = S(52);
    std::cout << std::fixed << std::setprecision(8) << ans << "\n";
    return 0;
}
