#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>
#include <future>

using namespace std;

const long long TARGET_MOD = 1000000007;

long long power(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128)res * base % mod;
        base = (__int128)base * base % mod;
        exp /= 2;
    }
    return res;
}

long long modInverse(long long n, long long mod) {
    return power(n, mod - 2, mod);
}

struct NTT_Mod {
    long long mod;
    long long root;
    long long root_pw;
};

const NTT_Mod P1 = {998244353, 3, 1 << 23};
const NTT_Mod P2 = {1004535809, 3, 1 << 21};
const NTT_Mod P3 = {469762049, 3, 1 << 26};

void ntt(vector<long long>& a, bool invert, const NTT_Mod& p) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        long long wlen = power(p.root, (p.mod - 1) / len, p.mod);
        if (invert) wlen = modInverse(wlen, p.mod);
        for (int i = 0; i < n; i += len) {
            long long w = 1;
            for (int j = 0; j < len / 2; j++) {
                long long u = a[i + j];
                long long v = (__int128)a[i + j + len / 2] * w % p.mod;
                a[i + j] = (u + v < p.mod) ? u + v : u + v - p.mod;
                a[i + j + len / 2] = (u - v >= 0) ? u - v : u - v + p.mod;
                w = (__int128)w * wlen % p.mod;
            }
        }
    }

    if (invert) {
        long long n_inv = modInverse(n, p.mod);
        for (long long& x : a) x = (__int128)x * n_inv % p.mod;
    }
}

vector<long long> convolve(const vector<long long>& a, const vector<long long>& b, const NTT_Mod& p) {
    vector<long long> fa(a.begin(), a.end());
    vector<long long> fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    ntt(fa, false, p);
    ntt(fb, false, p);
    for (int i = 0; i < n; i++) fa[i] = (__int128)fa[i] * fb[i] % p.mod;
    ntt(fa, true, p);

    vector<long long> result = fa;
    return result;
}

long long crt_combine(long long a1, long long a2, long long a3) {
    long long m1 = P1.mod, m2 = P2.mod, m3 = P3.mod;
    long long M = m1 * m2;
    long long inv_m1_m2 = modInverse(m1, m2);
    long long x = (a1 + (__int128)(a2 - a1 + m2) % m2 * inv_m1_m2 % m2 * m1);
    long long M_mod_m3 = M % m3;
    long long inv_M_m3 = modInverse(M_mod_m3, m3);
    long long k = (__int128)(a3 - (x % m3) + m3) % m3 * inv_M_m3 % m3;
    __int128 res = x + (__int128)k * M;
    return (long long)(res % TARGET_MOD);
}

vector<long long> threaded_convolution_crt(const vector<long long>& A, const vector<long long>& B) {
    auto f1 = async(launch::async, convolve, cref(A), cref(B), cref(P1));
    auto f2 = async(launch::async, convolve, cref(A), cref(B), cref(P2));
    auto f3 = async(launch::async, convolve, cref(A), cref(B), cref(P3));

    vector<long long> r1 = f1.get();
    vector<long long> r2 = f2.get();
    vector<long long> r3 = f3.get();

    int n = r1.size();
    vector<long long> res(n);

    int num_threads = thread::hardware_concurrency();
    vector<thread> threads;
    int chunk = n / num_threads + 1;
    
    for(int t=0; t<num_threads; ++t) {
        threads.emplace_back([&, t, chunk]() {
            int start = t * chunk;
            int end = min(start + chunk, n);
            for(int i=start; i<end; ++i) {
                res[i] = crt_combine(r1[i], r2[i], r3[i]);
            }
        });
    }
    for(auto& th : threads) th.join();
    
    return res;
}

class Solver {
    vector<long long> fact, invFact;
    int max_k;

public:
    Solver(int k) : max_k(k) {
        fact.resize(k + 1);
        invFact.resize(k + 1);
        fact[0] = 1;
        for (int i = 1; i <= k; i++) fact[i] = fact[i - 1] * i % TARGET_MOD;
        invFact[k] = modInverse(fact[k], TARGET_MOD);
        for (int i = k - 1; i >= 0; i--) invFact[i] = invFact[i + 1] * (i + 1) % TARGET_MOD;
    }

    long long nCr_pow16(int n, int r) {
        if (r < 0 || r > n) return 0;
        long long num = fact[n];
        long long den = (invFact[r] * invFact[n - r]) % TARGET_MOD;
        long long comb = (num * den) % TARGET_MOD;
        return power(comb, 16, TARGET_MOD); 
    }

    long long solve(int k, long long n) {
        long long phi = TARGET_MOD - 1;
        long long m_mod_phi = power(10, 11, phi); 
        long long m = n / k;
        long long m_exp = m % phi;
        
        long long lambda = (power(2, m_exp, TARGET_MOD) - 2 + TARGET_MOD) % TARGET_MOD;

        vector<long long> H(k + 1);
        for(int i=0; i<=k; ++i) {
            H[i] = power(invFact[i], 16, TARGET_MOD);
        }

        vector<long long> H2 = threaded_convolution_crt(H, H);

        vector<long long> U(k + 1);
        for(int s=0; s<=k; ++s) {
            long long A_s = (power(fact[s], 16, TARGET_MOD) * H2[s]) % TARGET_MOD;
            U[s] = (A_s * invFact[s]) % TARGET_MOD;
        }

        vector<long long> V(k + 1);
        long long lam_pow = 1;
        for(int j=0; j<=k; ++j) {
            V[j] = (lam_pow * invFact[j]) % TARGET_MOD;
            lam_pow = (lam_pow * lambda) % TARGET_MOD;
        }

        vector<long long> W = threaded_convolution_crt(U, V);

        long long ans = (fact[k] * W[k]) % TARGET_MOD;
        return ans;
    }
};

int main() {
    cout << "Running validation checks..." << endl;
    
    {
        Solver s(2);
        long long res = s.solve(2, 4);
        cout << "B(2, 4) = " << res << (res == 65550 ? " [PASS]" : " [FAIL]") << endl;
    }

    {
        Solver s(3);
        long long res = s.solve(3, 9);
        cout << "B(3, 9) = " << res << (res == 87273560 ? " [PASS]" : " [FAIL]") << endl;
    }

    cout << "-----------------------------------" << endl;
    
    int k_target = 100000;
    long long n_target = 10000000000000000LL;

    cout << "Solving for B(" << k_target << ", 10^16)..." << endl;
    
    Solver finalSolver(k_target); 
    
    auto start = chrono::high_resolution_clock::now();
    
    long long answer = finalSolver.solve(k_target, n_target);
    
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;

    cout << "Computation finished in " << diff.count() << " s" << endl;
    cout << "Answer: " << answer << endl;

    return 0;
}
