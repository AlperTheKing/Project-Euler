#include <iostream>
#include <vector>

using namespace std;

constexpr long long MOD = 1'000'000'007LL;
constexpr int kMaxBlock = 5;

long long mod_mul(long long a, long long b) {
    return static_cast<long long>((__int128)a * b % MOD);
}

long long power(long long base, long long exp) {
    long long res = 1 % MOD;
    base %= MOD;
    while (exp > 0) {
        if (exp & 1) res = mod_mul(res, base);
        base = mod_mul(base, base);
        exp >>= 1;
    }
    return res;
}

long long mod_inverse(long long n) {
    return power(n, MOD - 2);
}

void solve() {
    const int target = 10'000'000;

    const long long A[] = {0, 1, 2, 6, 18, 12};

    long long fact_small[kMaxBlock + 1];
    long long inv_fact_small[kMaxBlock + 1];
    long long coeff[kMaxBlock + 1];

    fact_small[0] = 1;
    inv_fact_small[0] = 1;
    coeff[0] = 0;
    for (int i = 1; i <= kMaxBlock; ++i) {
        fact_small[i] = mod_mul(fact_small[i - 1], i);
        inv_fact_small[i] = mod_inverse(fact_small[i]);
        coeff[i] = mod_mul(A[i], inv_fact_small[i]);
    }

    vector<long long> H(target + 1, 0);
    H[0] = 1;

    long long fact_n = 1;

    for (int n = 1; n <= target; ++n) {
        long long sum = 0;
        int limit = (n < kMaxBlock) ? n : kMaxBlock;
        for (int j = 1; j <= limit; ++j) {
            sum += mod_mul(coeff[j], H[n - j]);
            if (sum >= MOD) sum -= MOD;
        }
        H[n] = sum;
        fact_n = mod_mul(fact_n, n);
    }

    long long result = mod_mul(H[target], fact_n);
    cout << result << '\n';
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
