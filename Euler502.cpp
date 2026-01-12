#include <iostream>
#include <vector>
#include <algorithm>
#include <thread>
#include <future>
#include <chrono>

using namespace std;

long long MOD = 1000000007;

// --- Modular Arithmetic Utils ---
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

// --- Polynomial Class ---
struct Poly {
    vector<long long> coeffs;

    Poly() {}
    Poly(const vector<long long>& c) : coeffs(c) {}
    Poly(long long v) { if(v) coeffs = {v}; }

    void trim() {
        while (coeffs.size() > 0 && coeffs.back() == 0) {
            coeffs.pop_back();
        }
    }

    static Poly add(const Poly& A, const Poly& B) {
        Poly R;
        size_t n = max(A.coeffs.size(), B.coeffs.size());
        R.coeffs.resize(n, 0);
        for (size_t i = 0; i < n; i++) {
            long long a = (i < A.coeffs.size()) ? A.coeffs[i] : 0;
            long long b = (i < B.coeffs.size()) ? B.coeffs[i] : 0;
            long long sum = a + b;
            if (sum >= MOD) sum -= MOD;
            R.coeffs[i] = sum;
        }
        R.trim();
        return R;
    }

    static Poly sub(const Poly& A, const Poly& B) {
        Poly R;
        size_t n = max(A.coeffs.size(), B.coeffs.size());
        R.coeffs.resize(n, 0);
        for (size_t i = 0; i < n; i++) {
            long long a = (i < A.coeffs.size()) ? A.coeffs[i] : 0;
            long long b = (i < B.coeffs.size()) ? B.coeffs[i] : 0;
            long long diff = a - b;
            if (diff < 0) diff += MOD;
            R.coeffs[i] = diff;
        }
        R.trim();
        return R;
    }

    // Parallel Multiplication
    static Poly mul(const Poly& A, const Poly& B, int limit_deg = -1) {
        if (A.coeffs.empty() || B.coeffs.empty()) return Poly();
        size_t n = A.coeffs.size();
        size_t m = B.coeffs.size();
        
        size_t res_len = n + m - 1;
        if (limit_deg != -1) res_len = min(res_len, (size_t)limit_deg + 1);

        // Threshold for single thread
        if (n * m < 50000) {
            vector<long long> res(res_len, 0);
            for (size_t i = 0; i < n; i++) {
                if (A.coeffs[i] == 0) continue;
                for (size_t j = 0; j < m; j++) {
                    if (i + j >= res_len) break;
                    res[i + j] = (res[i + j] + A.coeffs[i] * B.coeffs[j]) % MOD;
                }
            }
            Poly R(res); R.trim(); return R;
        }

        int num_threads = thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 4;
        
        vector<vector<long long>> thread_res(num_threads, vector<long long>(res_len, 0));
        vector<thread> threads;

        for (int t = 0; t < num_threads; ++t) {
            threads.emplace_back([&, t]() {
                for (size_t i = t; i < n; i += num_threads) {
                    if (A.coeffs[i] == 0) continue;
                    for (size_t j = 0; j < m; j++) {
                        if (i + j >= res_len) break;
                        thread_res[t][i + j] = (thread_res[t][i + j] + A.coeffs[i] * B.coeffs[j]) % MOD;
                    }
                }
            });
        }

        for (auto& th : threads) th.join();

        Poly R;
        R.coeffs.resize(res_len, 0);
        for (size_t i = 0; i < res_len; ++i) {
            long long sum = 0;
            for (int t = 0; t < num_threads; ++t) {
                sum += thread_res[t][i];
            }
            R.coeffs[i] = sum % MOD;
        }
        R.trim();
        return R;
    }
};

// --- Matrix Class ---
struct Matrix {
    Poly mat[2][2];
    int trunc_deg; // -1 for exact
    
    Matrix(int deg = -1) : trunc_deg(deg) {}

    static Matrix mul(const Matrix& A, const Matrix& B) {
        Matrix C(A.trunc_deg);
        for(int i=0; i<2; i++) {
            for(int k=0; k<2; k++) {
                if (A.mat[i][k].coeffs.empty()) continue;
                for(int j=0; j<2; j++) {
                    if (B.mat[k][j].coeffs.empty()) continue;
                    Poly prod = Poly::mul(A.mat[i][k], B.mat[k][j], A.trunc_deg);
                    C.mat[i][j] = Poly::add(C.mat[i][j], prod);
                }
            }
        }
        return C;
    }

    static Matrix identity(int deg = -1) {
        Matrix I(deg);
        I.mat[0][0] = Poly(1);
        I.mat[1][1] = Poly(1);
        return I;
    }
    
    static Matrix power(Matrix base, long long exp) {
        Matrix res = Matrix::identity(base.trunc_deg);
        while (exp > 0) {
            if (exp & 1) res = Matrix::mul(res, base);
            base = Matrix::mul(base, base);
            exp >>= 1;
        }
        return res;
    }
};

// --- Core Logic ---

// Get transition matrix M
// type=1 (S):  [[1+x, x], [-x, 1-x]]
// type=-1 (D): [[-(1+x), -x], [-x, 1-x]]
Matrix get_matrix(int type, int trunc_deg) {
    Matrix M(trunc_deg);
    Poly One({1});
    Poly X({0, 1});
    Poly OnePlusX({1, 1});
    Poly OneMinusX({1, MOD-1});
    Poly NegX({0, MOD-1});
    
    if (type == 1) {
        M.mat[0][0] = OnePlusX;
        M.mat[0][1] = X;
        M.mat[1][0] = NegX;
        M.mat[1][1] = OneMinusX;
    } else {
        M.mat[0][0] = Poly::sub(Poly(), OnePlusX); // -(1+x)
        M.mat[0][1] = NegX;
        M.mat[1][0] = NegX;
        M.mat[1][1] = OneMinusX;
    }
    return M;
}

// Solver 1: Bostan-Mori for large W, small H
long long solve_bostan_mori(long long w, long long h, int type) {
    Matrix M = get_matrix(type, -1);
    Matrix Mh = Matrix::power(M, h);
    
    // Result vector = Mh * [0, 1]^T
    Poly P = Mh.mat[0][1];
    Poly Q = Mh.mat[1][1];
    
    // Extract [x^w] P/Q
    long long k = w;
    while (k > 0) {
        Poly Q_neg = Q;
        for (size_t i = 1; i < Q_neg.coeffs.size(); i += 2) {
            if (Q_neg.coeffs[i] != 0) Q_neg.coeffs[i] = MOD - Q_neg.coeffs[i];
        }

        Poly U = Poly::mul(P, Q_neg);
        Poly V = Poly::mul(Q, Q_neg);

        vector<long long> next_P, next_Q;
        int parity = (k % 2);
        for (size_t i = parity; i < U.coeffs.size(); i += 2) next_P.push_back(U.coeffs[i]);
        for (size_t i = 0; i < V.coeffs.size(); i += 2) next_Q.push_back(V.coeffs[i]);

        P = Poly(next_P);
        Q = Poly(next_Q);
        k /= 2;
    }
    
    long long p0 = (P.coeffs.empty()) ? 0 : P.coeffs[0];
    long long q0 = (Q.coeffs.empty()) ? 0 : Q.coeffs[0];
    return p0 * modInverse(q0) % MOD;
}

// Solver 2: Matrix Exponentiation with truncation for large H, small W
long long solve_trunc(long long w, long long h, int type) {
    Matrix M = get_matrix(type, w);
    Matrix Mh = Matrix::power(M, h);
    
    Poly P = Mh.mat[0][1];
    Poly Q = Mh.mat[1][1];
    
    // Compute P * Q^-1 mod x^{w+1}
    // Q[0] is 1 because det(M)=1/-1 and traces -> (1-x) dominant term
    // Let's verify Q[0] != 0.
    if (Q.coeffs.empty() || Q.coeffs[0] == 0) return 0;
    
    vector<long long> invQ(w + 1, 0);
    long long invQ0 = modInverse(Q.coeffs[0]);
    invQ[0] = invQ0;
    
    for (int i = 1; i <= w; ++i) {
        long long sum = 0;
        for (int j = 1; j <= i; ++j) {
            if (j < Q.coeffs.size()) {
                sum = (sum + Q.coeffs[j] * invQ[i - j]) % MOD;
            }
        }
        invQ[i] = (MOD - sum) * invQ0 % MOD;
    }
    
    long long ans = 0;
    for (int i = 0; i <= w; ++i) {
        long long p_val = (i < P.coeffs.size()) ? P.coeffs[i] : 0;
        if (p_val != 0) {
            ans = (ans + p_val * invQ[w - i]) % MOD;
        }
    }
    return ans;
}

long long get_count_le(long long w, long long h, int type) {
    if (h == 0) return 0;
    // Heuristic switch
    if (w > 20000) {
        return solve_bostan_mori(w, h, type);
    } else {
        return solve_trunc(w, h, type);
    }
}

// Calculate Even Count for height exactly H
// Even(h) = (S_h + D_h)/2
// Result = Even(h) - Even(h-1)
long long solve_exact_h(long long w, long long h) {
    if (h <= 0) return 0;
    long long S_h = get_count_le(w, h, 1);
    long long D_h = get_count_le(w, h, -1);
    long long even_h = (S_h + D_h) % MOD * modInverse(2) % MOD;
    
    long long S_h_1 = get_count_le(w, h - 1, 1);
    long long D_h_1 = get_count_le(w, h - 1, -1);
    long long even_h_1 = (S_h_1 + D_h_1) % MOD * modInverse(2) % MOD;
    
    return (even_h - even_h_1 + MOD) % MOD;
}

int main() {
    // --- Validation ---
    cout << "Validating..." << endl;
    
    long long v1 = solve_exact_h(4, 2);
    cout << "F(4, 2) = " << v1 << " (Expected 10)" << endl;
    
    long long v2 = solve_exact_h(13, 10);
    cout << "F(13, 10) = " << v2 << " (Expected " << 3729050610636LL % MOD << ")" << endl;
    
    long long v3 = solve_exact_h(10, 13);
    cout << "F(10, 13) = " << v3 << " (Expected " << 37959702514LL % MOD << ")" << endl;
    
    if (v1 != 10 || v2 != 3729050610636LL % MOD || v3 != 37959702514LL % MOD) {
        cout << "VALIDATION FAILED" << endl;
        return 1;
    }
    cout << "Validation Passed." << endl;
    cout << "------------------" << endl;

    // --- Solution ---
    long long ans1 = solve_exact_h(1000000000000LL, 100);
    cout << "Term 1: " << ans1 << endl;
    
    long long ans2 = solve_exact_h(10000, 10000);
    cout << "Term 2: " << ans2 << endl;
    
    long long ans3 = solve_exact_h(100, 1000000000000LL);
    cout << "Term 3: " << ans3 << endl;

    long long total = (ans1 + ans2 + ans3) % MOD;
    cout << "Final Answer: " << total << endl;

    return 0;
}