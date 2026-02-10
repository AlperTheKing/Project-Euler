#include <cstdint>
#include <iostream>
#include <vector>

// Project Euler 602: Expected product of head-counts
//
// Let n be the number of friends (excluding Alice). The game ends when Alice gets her
// first Head. Let r be the number of Alice's Tails before that Head. Then:
//   P(r = t) = p^t (1-p),  t >= 0,
// and each friend makes exactly r tosses, with head probability (1-p).
//
// For a fixed r, each friend's head count H ~ Binomial(r, 1-p), so E[H | r] = r(1-p).
// The friends are independent given r, hence:
//   e(n,p) = E[ prod_i H_i ] = E[ (E[H | r])^n ] = sum_{r>=0} p^r (1-p) * (r(1-p))^n
//          = (1-p)^{n+1} * sum_{r>=0} r^n p^r.
//
// The known generating function identity is:
//   sum_{r>=0} r^n x^r = x * A_n(x) / (1-x)^{n+1},
// where A_n(x) is the Eulerian polynomial:
//   A_n(x) = sum_{m=0..n-1} <n,m> x^m,
// and <n,m> is the Eulerian number (permutations of {1..n} with m descents).
//
// Therefore:
//   e(n,p) = (1-p)^{n+1} * p * A_n(p) / (1-p)^{n+1} = p * A_n(p),
// so the coefficient c(n,k) of p^k in e(n,p) is:
//   c(n,k) = <n, k-1>.
//
// We compute Eulerian(n,m) modulo MOD using:
//   <n,m> = sum_{j=0..m+1} (-1)^j * C(n+1, j) * (m+1-j)^n.

using u64 = std::uint64_t;

static constexpr int MOD = 1000000007;

static inline int addmod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}
static inline int submod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}
static inline int mulmod(long long a, long long b) {
    return (int)((a * b) % MOD);
}

static int mod_pow(int a, u64 e) {
    long long r = 1;
    long long x = a;
    while (e) {
        if (e & 1) r = (r * x) % MOD;
        x = (x * x) % MOD;
        e >>= 1;
    }
    return (int)r;
}

static int eulerian_mod(u64 n, u64 m) {
    // Requires 0 <= m <= n-1, n >= 1.
    const u64 N = n + 1; // binomial top
    const u64 JMAX = m + 1;

    std::vector<int> inv(JMAX + 2);
    if (JMAX + 1 >= 1) inv[1] = 1;
    for (u64 i = 2; i <= JMAX + 1; ++i) {
        inv[i] = (int)(MOD - (long long)(MOD / (int)i) * inv[MOD % (int)i] % MOD);
    }

    int ans = 0;
    long long binom = 1; // C(N,0)
    for (u64 j = 0; j <= JMAX; ++j) {
        const int base = (int)(JMAX - j); // (m+1-j)
        const int pown = mod_pow(base, n);
        const int term = mulmod((int)binom, pown);
        ans = (j & 1) ? submod(ans, term) : addmod(ans, term);

        if (j == JMAX) break;
        // C(N, j+1) = C(N,j) * (N-j)/(j+1).
        binom = binom * (long long)(N - j) % MOD;
        binom = binom * inv[j + 1] % MOD;
    }
    return ans;
}

static int c_mod(u64 n, u64 k) {
    // coefficient of p^k in e(n,p) is Eulerian(n,k-1)
    if (k == 0 || k > n) return 0;
    return eulerian_mod(n, k - 1);
}

int main() {
    // Statement examples/validations.
    if (c_mod(3, 1) != 1 || c_mod(3, 2) != 4 || c_mod(3, 3) != 1) {
        std::cerr << "Validation failed: e(3,p)\n";
        return 1;
    }
    if (c_mod(100, 40) != 986699437) {
        std::cerr << "Validation failed: c(100,40)\n";
        return 1;
    }

    std::cout << c_mod(10000000ULL, 4000000ULL) << "\n";
    return 0;
}

