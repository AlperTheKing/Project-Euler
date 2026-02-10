#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = long long;

static constexpr int MOD = 1'008'691'207;

static inline int mod_norm(i64 x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return (int)x;
}

static int f_formula(i64 n) {
    if (n <= 2) return 0;

    // f(n)=n!-2(n-1)!-(n-2)!-(n-3)!+2+(n-3)*sum_{k=0}^{n-4} k! (mod MOD).
    i64 fact = 1;
    i64 sum_fact = (n >= 4 ? 1 : 0);  // 0!

    const i64 upto = n - 4;
    for (i64 i = 1; i <= upto; ++i) {
        fact = (fact * i) % MOD;
        sum_fact += fact;
        if (sum_fact >= MOD) sum_fact -= MOD;
    }

    i64 fn3 = (n == 3 ? 1 : 0), fn2 = 0, fn1 = 0, fn = 0;
    for (i64 i = std::max<i64>(1, upto + 1); i <= n; ++i) {
        fact = (fact * i) % MOD;
        if (i == n - 3) fn3 = fact;
        if (i == n - 2) fn2 = fact;
        if (i == n - 1) fn1 = fact;
        if (i == n) fn = fact;
    }

    i64 res = fn;
    res = (res - 2 * fn1) % MOD;
    res = (res - fn2) % MOD;
    res = (res - fn3) % MOD;
    res = (res + 2) % MOD;
    res = (res + ((n - 3) % MOD) * sum_fact) % MOD;
    return mod_norm(res);
}

static bool is_open(const std::vector<int> &p) {
    const int n = (int)p.size();
    std::vector<std::vector<char>> blocked(n, std::vector<char>(n, 0));
    for (int r = 0; r < n; ++r) blocked[r][p[r]] = 1;
    if (blocked[0][0] || blocked[n - 1][n - 1]) return false;

    std::vector<std::vector<char>> dp(n, std::vector<char>(n, 0));
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (blocked[r][c]) continue;
            if (r == 0 && c == 0) {
                dp[r][c] = 1;
            } else {
                dp[r][c] = (char)((r > 0 && dp[r - 1][c]) || (c > 0 && dp[r][c - 1]));
            }
        }
    }
    return dp[n - 1][n - 1];
}

static i64 brute_count(int n) {
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    i64 cnt = 0;
    do {
        if (is_open(p)) ++cnt;
    } while (std::next_permutation(p.begin(), p.end()));
    return cnt;
}

int main() {
    assert(f_formula(3) == 2);
    assert(f_formula(5) == 70);
    for (int n = 3; n <= 9; ++n) assert((i64)f_formula(n) == brute_count(n));

    std::cout << f_formula(100'000'000) << "\n";
    return 0;
}

