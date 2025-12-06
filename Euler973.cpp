#include <iostream>
#include <vector>
#include <future>
#include <cstdlib>

using namespace std;

long long MOD = 1e9 + 7;

long long power(long long b, long long e) {
    long long r = 1;
    b %= MOD;
    while (e) {
        if (e & 1) r = (r * b) % MOD;
        b = (b * b) % MOD;
        e >>= 1;
    }
    return r;
}

long long get_stable(int k, int n) {
    long long s = 1LL << k;
    if (n == s) return 1;
    long long t1 = (n - s + 2) % MOD;
    long long exp = n - s - 1;
    return (t1 * power(2, exp)) % MOD;
}

long long solve_bit(int k, int n) {
    long long s = 1LL << k;
    if (n < s) return 0;
    long long es = 1LL << (k + 1);
    if (n < es) return get_stable(k, n);
    vector<long long> dp(n + 1, 0);
    for (int i = s; i < es; ++i) dp[i] = get_stable(k, i);
    long long c_sum = 0;
    long long ub = 1LL << k;
    if (ub > 2) {
        for (int i = es - 2; i >= es - ub + 1; --i)
            c_sum = (c_sum + dp[i]) % MOD;
    }
    for (int i = es; i <= n; ++i) {
        long long t1 = (3 * dp[i - 1]) % MOD;
        long long t3 = (4 * dp[i - ub]) % MOD;
        long long t4 = (4 * dp[i - ub - 1]) % MOD;
        long long val = (t1 - c_sum - t3 + t4) % MOD;
        dp[i] = (val + 2 * MOD) % MOD;
        if (ub > 2) {
            c_sum = (c_sum + dp[i - 1]) % MOD;
            c_sum = (c_sum - dp[i - ub + 1] + MOD) % MOD;
        }
    }
    return dp[n];
}

long long solve_X(int n) {
    long long tot = 0;
    if (n % 2 != 0) tot = (power(2, n - 1) - 1 + MOD) % MOD;
    vector<future<long long>> futs;
    for (int k = 1; k < 20; ++k) {
        if ((1LL << k) > n) break;
        futs.push_back(async(launch::async, solve_bit, k, n));
    }
    for (int k = 0; k < futs.size(); ++k) {
        long long val = futs[k].get();
        long long term = (val * power(2, k + 1)) % MOD;
        tot = (tot + term) % MOD;
    }
    return tot;
}

int main() {
    if (solve_X(2) != 2) exit(1);
    if (solve_X(4) != 14) exit(2);
    if (solve_X(10) != 1418) exit(3);
    cout << solve_X(10000) << endl;
    return 0;
}