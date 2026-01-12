#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

using namespace std;

struct Matrix {
    int n = 0;
    vector<long double> a;

    Matrix() = default;
    explicit Matrix(int n_) : n(n_), a(static_cast<size_t>(n_) * n_, 0.0L) {}

    long double& operator()(int r, int c) {
        return a[static_cast<size_t>(r) * n + c];
    }

    long double operator()(int r, int c) const {
        return a[static_cast<size_t>(r) * n + c];
    }
};

static Matrix identity_matrix(int n) {
    Matrix I(n);
    for (int i = 0; i < n; ++i) {
        I(i, i) = 1.0L;
    }
    return I;
}

static Matrix multiply(const Matrix& A, const Matrix& B) {
    int n = A.n;
    Matrix C(n);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            long double aik = A(i, k);
            if (aik == 0.0L) {
                continue;
            }
            const long double* brow = &B.a[static_cast<size_t>(k) * n];
            long double* crow = &C.a[static_cast<size_t>(i) * n];
            for (int j = 0; j < n; ++j) {
                crow[j] += aik * brow[j];
            }
        }
    }
    return C;
}

static vector<long double> multiply_vec(const Matrix& A, const vector<long double>& v) {
    int n = A.n;
    vector<long double> res(n, 0.0L);
    for (int i = 0; i < n; ++i) {
        long double sum = 0.0L;
        const long double* row = &A.a[static_cast<size_t>(i) * n];
        for (int j = 0; j < n; ++j) {
            sum += row[j] * v[j];
        }
        res[i] = sum;
    }
    return res;
}

static vector<Matrix> build_pow2(const Matrix& R, int bits) {
    vector<Matrix> pow2;
    pow2.reserve(bits);
    pow2.push_back(R);
    for (int i = 1; i < bits; ++i) {
        pow2.push_back(multiply(pow2.back(), pow2.back()));
    }
    return pow2;
}

static Matrix pow_from_bits(const vector<Matrix>& pow2, const vector<int>& bits, int n) {
    Matrix res = identity_matrix(n);
    for (int bit : bits) {
        res = multiply(res, pow2[bit]);
    }
    return res;
}

static bool solve_linear_system(Matrix A, vector<long double> b, vector<long double>& x) {
    int n = A.n;
    x.assign(n, 0.0L);
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        long double max_abs = fabsl(A(i, i));
        for (int r = i + 1; r < n; ++r) {
            long double val = fabsl(A(r, i));
            if (val > max_abs) {
                max_abs = val;
                pivot = r;
            }
        }
        if (max_abs < 1e-30L) {
            return false;
        }
        if (pivot != i) {
            for (int c = i; c < n; ++c) {
                swap(A(i, c), A(pivot, c));
            }
            swap(b[i], b[pivot]);
        }
        long double diag = A(i, i);
        for (int c = i; c < n; ++c) {
            A(i, c) /= diag;
        }
        b[i] /= diag;
        for (int r = 0; r < n; ++r) {
            if (r == i) {
                continue;
            }
            long double factor = A(r, i);
            if (factor == 0.0L) {
                continue;
            }
            for (int c = i; c < n; ++c) {
                A(r, c) -= factor * A(i, c);
            }
            b[r] -= factor * b[i];
        }
    }
    x = b;
    return true;
}

struct Transition {
    int from = 0;
    int to = 0;
    int exp_index = 0;
    long double prob = 0.0L;
};

struct ExpInfo {
    long long exp = 0;
    vector<int> bits;
};

class LotterySolver {
public:
    LotterySolver(long long m_val, int max_n)
        : m(m_val), d(static_cast<int>(m_val - 1)), N(max_n) {
        build_transitions();
    }

    long double probability_no_ruin(long long s) {
        if (s < m) {
            return 0.0L;
        }
        Matrix R = compute_R();
        vector<long double> y0 = compute_y0(R);
        long long c = s - d;
        long long k = (c - 1) / d;
        int r = static_cast<int>((c - 1) % d);
        if (k == 0) {
            return 1.0L - y0[r];
        }
        Matrix Rk = matrix_power(R, k);
        vector<long double> yk = multiply_vec(Rk, y0);
        return 1.0L - yk[r];
    }

private:
    long long m;
    int d;
    int N;

    vector<long long> pow2n;
    vector<long double> probs;
    Matrix A_minus;
    vector<Transition> transitions;
    vector<ExpInfo> exp_info;
    void build_transitions() {
        pow2n.resize(N);
        probs.resize(N);
        for (int n = 0; n < N; ++n) {
            pow2n[n] = 1LL << n;
            probs[n] = ldexp(1.0L, -(n + 1));
        }

        A_minus = Matrix(d);
        vector<long long> exps;
        transitions.clear();

        for (int n = 0; n < N; ++n) {
            long long delta = pow2n[n] - m;
            long double prob = probs[n];
            for (int r = 0; r < d; ++r) {
                long long c = d + r + 1;
                long long cprime = c + delta;
                long long kprime = (cprime - 1) / d;
                long long i = kprime - 1;
                int rprime = static_cast<int>((cprime - 1) % d);
                if (i == -1) {
                    A_minus(r, rprime) += prob;
                } else {
                    long long exp = i + 1;
                    exps.push_back(exp);
                    transitions.push_back({r, rprime, 0, prob});
                }
            }
        }

        sort(exps.begin(), exps.end());
        exps.erase(unique(exps.begin(), exps.end()), exps.end());

        unordered_map<long long, int> exp_index;
        exp_info.clear();
        exp_info.reserve(exps.size());
        for (size_t i = 0; i < exps.size(); ++i) {
            exp_index[exps[i]] = static_cast<int>(i);
            ExpInfo info;
            info.exp = exps[i];
            long long e = exps[i];
            int bit = 0;
            while (e > 0) {
                if (e & 1LL) {
                    info.bits.push_back(bit);
                }
                e >>= 1LL;
                ++bit;
            }
            exp_info.push_back(info);
        }

        size_t idx = 0;
        for (int n = 0; n < N; ++n) {
            long long delta = pow2n[n] - m;
            long double prob = probs[n];
            for (int r = 0; r < d; ++r) {
                long long c = d + r + 1;
                long long cprime = c + delta;
                long long kprime = (cprime - 1) / d;
                long long i = kprime - 1;
                if (i == -1) {
                    continue;
                }
                long long exp = i + 1;
                transitions[idx].exp_index = exp_index[exp];
                transitions[idx].prob = prob;
                ++idx;
            }
        }
    }

    Matrix matrix_power(const Matrix& R, long long exp) const {
        if (exp == 0) {
            return identity_matrix(d);
        }
        int max_bit = 0;
        long long temp = exp;
        while (temp > 0) {
            ++max_bit;
            temp >>= 1LL;
        }
        vector<Matrix> pow2 = build_pow2(R, max_bit);
        vector<int> bits;
        long long e = exp;
        int bit = 0;
        while (e > 0) {
            if (e & 1LL) {
                bits.push_back(bit);
            }
            e >>= 1LL;
            ++bit;
        }
        return pow_from_bits(pow2, bits, d);
    }

    Matrix compute_R() {
        Matrix R = A_minus;
        const long double tol = 1e-18L;
        const int max_iter = 3000;

        int max_bit = 0;
        for (const auto& info : exp_info) {
            if (!info.bits.empty()) {
                max_bit = max(max_bit, info.bits.back() + 1);
            }
        }
        if (max_bit == 0) {
            max_bit = 1;
        }

        for (int iter = 0; iter < max_iter; ++iter) {
            vector<Matrix> pow2 = build_pow2(R, max_bit);
            vector<Matrix> Rexp(exp_info.size(), Matrix(d));
            for (size_t i = 0; i < exp_info.size(); ++i) {
                Rexp[i] = pow_from_bits(pow2, exp_info[i].bits, d);
            }

            Matrix Rnew = A_minus;
            for (const auto& tr : transitions) {
                const Matrix& mat = Rexp[tr.exp_index];
                const long double* src_row = &mat.a[static_cast<size_t>(tr.to) * d];
                long double* dst_row = &Rnew.a[static_cast<size_t>(tr.from) * d];
                for (int j = 0; j < d; ++j) {
                    dst_row[j] += tr.prob * src_row[j];
                }
            }

            long double diff = 0.0L;
            for (size_t i = 0; i < R.a.size(); ++i) {
                diff = max(diff, fabsl(Rnew.a[i] - R.a[i]));
            }
            R = Rnew;
            if (diff < tol) {
                break;
            }
            if (iter == max_iter - 1) {
                cerr << "Warning: R iteration did not fully converge (diff=" << setprecision(6) << diff << ")\n";
            }
        }
        return R;
    }

    vector<long double> compute_y0(const Matrix& R) {
        vector<long long> boundary_exps;
        boundary_exps.reserve(static_cast<size_t>(2 * N));

        for (int n = 0; n < N; ++n) {
            long long delta = pow2n[n] - m;
            for (int r = 0; r < d; ++r) {
                long long c = r + 1;
                long long cprime = c + delta;
                if (cprime <= 0) {
                    continue;
                }
                long long kprime = (cprime - 1) / d;
                if (kprime >= 1) {
                    boundary_exps.push_back(kprime);
                }
            }
        }

        sort(boundary_exps.begin(), boundary_exps.end());
        boundary_exps.erase(unique(boundary_exps.begin(), boundary_exps.end()), boundary_exps.end());

        unordered_map<long long, int> exp_index;
        vector<vector<int>> exp_bits(boundary_exps.size());
        int max_bit = 0;
        for (size_t i = 0; i < boundary_exps.size(); ++i) {
            exp_index[boundary_exps[i]] = static_cast<int>(i);
            long long e = boundary_exps[i];
            int bit = 0;
            while (e > 0) {
                if (e & 1LL) {
                    exp_bits[i].push_back(bit);
                }
                e >>= 1LL;
                ++bit;
            }
            if (!exp_bits[i].empty()) {
                max_bit = max(max_bit, exp_bits[i].back() + 1);
            }
        }
        if (max_bit == 0) {
            max_bit = 1;
        }

        vector<Matrix> pow2 = build_pow2(R, max_bit);
        vector<Matrix> Rexp(boundary_exps.size(), Matrix(d));
        for (size_t i = 0; i < boundary_exps.size(); ++i) {
            Rexp[i] = pow_from_bits(pow2, exp_bits[i], d);
        }

        Matrix M(d);
        vector<long double> b(d, 0.0L);

        for (int n = 0; n < N; ++n) {
            long long delta = pow2n[n] - m;
            long double prob = probs[n];
            for (int r = 0; r < d; ++r) {
                long long c = r + 1;
                long long cprime = c + delta;
                if (cprime <= 0) {
                    b[r] += prob;
                    continue;
                }
                long long kprime = (cprime - 1) / d;
                int rprime = static_cast<int>((cprime - 1) % d);
                if (kprime == 0) {
                    M(r, rprime) += prob;
                } else {
                    int idx = exp_index[kprime];
                    const Matrix& mat = Rexp[idx];
                    const long double* src_row = &mat.a[static_cast<size_t>(rprime) * d];
                    long double* dst_row = &M.a[static_cast<size_t>(r) * d];
                    for (int j = 0; j < d; ++j) {
                        dst_row[j] += prob * src_row[j];
                    }
                }
            }
        }

        Matrix A = identity_matrix(d);
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                A(i, j) -= M(i, j);
            }
        }

        vector<long double> y0;
        if (!solve_linear_system(A, b, y0)) {
            cerr << "Failed to solve boundary system.\n";
            y0.assign(d, 1.0L);
        }
        return y0;
    }
};

static void run_check(long long m, long long s, long double expected, long double tol) {
    const int N = 50;
    LotterySolver solver(m, N);
    long double res = solver.probability_no_ruin(s);
    long double diff = fabsl(res - expected);
    cout << "p_" << m << "(" << s << ") = " << fixed << setprecision(7) << res;
    if (diff <= tol) {
        cout << " [PASS]";
    } else {
        cout << " [FAIL]";
    }
    cout << "\n";
}

int main() {
    cout << "--- Validation Checkpoints ---\n";
    run_check(2, 2, 0.2522L, 5e-4L);
    run_check(2, 5, 0.6873L, 5e-4L);
    run_check(6, 10000, 0.9952L, 5e-4L);

    cout << "\n--- Final Solution ---\n";
    const long long m = 15;
    const long long s = 1000000000LL;
    const int N = 50;
    LotterySolver solver(m, N);
    long double ans = solver.probability_no_ruin(s);
    cout << fixed << setprecision(7) << ans << "\n";
    return 0;
}
