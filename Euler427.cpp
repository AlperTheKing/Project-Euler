#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

long long MOD = 1000000009;

long long power(long long base, long long exp) {
    long long res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp /= 2;
    }
    return res;
}

long long modInverse(long long n) {
    return power(n, MOD - 2);
}

vector<long long> fact;
vector<long long> invFact;
vector<long long> powN;
vector<long long> powN_1; 

void precompute(int n) {
    fact.resize(n + 1);
    invFact.resize(n + 1);
    powN.resize(n + 1);
    powN_1.resize(n + 1);

    fact[0] = 1;
    powN[0] = 1;
    powN_1[0] = 1;
    long long n_mod = n % MOD;
    long long n_1_mod = (MOD - (n - 1) % MOD) % MOD; 

    for (int i = 1; i <= n; i++) {
        fact[i] = (fact[i - 1] * i) % MOD;
        powN[i] = (powN[i - 1] * n_mod) % MOD;
        powN_1[i] = (powN_1[i - 1] * n_1_mod) % MOD;
    }

    invFact[n] = modInverse(fact[n]);
    for (int i = n - 1; i >= 0; i--) {
        invFact[i] = (invFact[i + 1] * (i + 1)) % MOD;
    }
}

long long nCr(int n, int r) {
    if (r < 0 || r > n) return 0;
    return fact[n] * invFact[r] % MOD * invFact[n - r] % MOD;
}

long long get_C(int M, int k, int n_val) {
    if (M < 0) return 0;
    long long sum = 0;
    int max_p = M / (k + 1);
    
    for (int p = 0; p <= max_p; ++p) {
        long long term = nCr(M - p * k, p);
        term = (term * powN[M - p * (k + 1)]) % MOD;
        term = (term * powN_1[p]) % MOD;
        sum = (sum + term) % MOD;
    }
    return sum;
}

long long solve_f(int n) {
    precompute(n);
    
    long long total_sequences = power(n, n + 1); 
    long long sum_Nk = 0;

    #pragma omp parallel for reduction(+:sum_Nk)
    for (int k = 1; k < n; ++k) {
        long long c_n_1 = get_C(n - 1, k, n);
        long long c_n_1_k = get_C(n - 1 - k, k, n);
        
        long long W_n = (n * c_n_1) % MOD;
        long long sub = (n * c_n_1_k) % MOD;
        
        W_n = (W_n - sub + MOD) % MOD;
        
        sum_Nk = (sum_Nk + W_n) % MOD;
    }

    long long ans = (total_sequences - sum_Nk + MOD) % MOD;
    return ans;
}

int main() {
    struct TestCase {
        int n;
        long long expected;
    };
    
    vector<TestCase> tests = {
        {3, 45},
        {7, 1403689},
        {11, 496890792} 
    };

    bool all_passed = true;
    for (const auto& t : tests) {
        long long result = solve_f(t.n);
        if (result != t.expected) {
            cout << "Test failed for n=" << t.n << ". Expected " << t.expected << ", got " << result << endl;
            all_passed = false;
        } else {
            cout << "Test passed for n=" << t.n << endl;
        }
    }

    if (all_passed) {
        int target_n = 7500000;
        cout << "Computing f(" << target_n << ")..." << endl;
        long long final_ans = solve_f(target_n);
        cout << "f(" << target_n << ") mod " << MOD << " = " << final_ans << endl;
    }

    return 0;
}