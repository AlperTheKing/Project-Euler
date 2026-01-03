#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

using namespace std;

namespace {

struct Precomp {
    int a = 0;
    int b = 0;
    int L = 0;
    int M = 0;
    vector<int> left_cnt;
    vector<int> right_cnt;
    vector<int> left_k1;
    vector<int> left_k2;
    vector<int> right_k1;
    vector<int> right_k2;
    vector<uint8_t> mid;
    vector<double> left_sum;
};

Precomp build_precomp(int n, int m) {
    Precomp p;
    p.a = n + 5;
    p.b = m + 5;
    p.L = m - n + 1;
    p.M = p.b;

    int size = p.M + 1;
    p.left_cnt.assign(size, 0);
    p.right_cnt.assign(size, 0);
    p.left_k1.assign(size, 1);
    p.left_k2.assign(size, 0);
    p.right_k1.assign(size, 1);
    p.right_k2.assign(size, 0);
    p.mid.assign(size, 0);
    p.left_sum.assign(size, 0.0);

    for (int r = 1; r <= p.M; ++r) {
        int left_l = p.a;
        int left_r = min(p.b, r - 1);
        if (left_r >= left_l) {
            int cnt = left_r - left_l + 1;
            p.left_cnt[r] = cnt;
            p.left_sum[r] = 0.5 * (left_l + left_r) * cnt;
            p.left_k1[r] = max(1, r - left_r);
            p.left_k2[r] = r - left_l;
        }

        int right_l = max(p.a, r + 1);
        int right_r = p.b;
        if (right_r >= right_l) {
            int cnt = right_r - right_l + 1;
            p.right_cnt[r] = cnt;
            p.right_k1[r] = max(1, right_l - r);
            p.right_k2[r] = right_r - r;
        }

        if (r >= p.a && r <= p.b) p.mid[r] = 1;
    }

    return p;
}

struct EResult {
    double value = 0.0;
    bool converged = false;
    int iterations = 0;
};

// Dense Gaussian elimination with partial pivoting.
bool solve_linear(vector<double>& mat, vector<double>& rhs, int n, vector<double>& sol) {
    sol.assign(n, 0.0);
    const double eps = 1e-14;

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        double maxabs = fabs(mat[i * n + i]);
        for (int row = i + 1; row < n; ++row) {
            double val = fabs(mat[row * n + i]);
            if (val > maxabs) {
                maxabs = val;
                pivot = row;
            }
        }
        if (!isfinite(maxabs) || maxabs < eps) return false;

        if (pivot != i) {
            for (int col = i; col < n; ++col) {
                swap(mat[i * n + col], mat[pivot * n + col]);
            }
            swap(rhs[i], rhs[pivot]);
        }

        double pivot_val = mat[i * n + i];
        for (int row = i + 1; row < n; ++row) {
            double factor = mat[row * n + i] / pivot_val;
            if (factor == 0.0) continue;
            for (int col = i; col < n; ++col) {
                mat[row * n + col] -= factor * mat[i * n + col];
            }
            rhs[row] -= factor * rhs[i];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = rhs[i];
        for (int col = i + 1; col < n; ++col) {
            sum -= mat[i * n + col] * sol[col];
        }
        sol[i] = sum / mat[i * n + i];
        if (!isfinite(sol[i])) return false;
    }
    return true;
}

EResult compute_E(int n, int m, double tol = 1e-12, int max_iter = 8000) {
    (void)tol;
    (void)max_iter;
    Precomp p = build_precomp(n, m);
    const int M = p.M;
    const int a = p.a;
    const int b = p.b;
    const int L = p.L;
    const double inv_L = 1.0 / static_cast<double>(L);
    const int N = 2 * M + 2;

    vector<double> mat(static_cast<size_t>(N) * N, 0.0);
    vector<double> rhs(N, 0.0);

    auto idx0 = [M](int r) { return r - 1; };
    auto idx1 = [M](int r) { return M + (r - 1); };
    const int idxAvg0 = 2 * M;
    const int idxAvg1 = 2 * M + 1;

    for (int r = 1; r <= M; ++r) {
        int row0 = idx0(r);
        int row1 = idx1(r);

        mat[row0 * N + row0] = 1.0;
        if (p.mid[r]) mat[row0 * N + idxAvg0] -= inv_L;
        if (p.left_cnt[r] > 0) {
            for (int k = p.left_k1[r]; k <= p.left_k2[r]; ++k) {
                mat[row0 * N + idx1(k)] -= inv_L;
            }
        }
        if (p.right_cnt[r] > 0) {
            for (int k = p.right_k1[r]; k <= p.right_k2[r]; ++k) {
                mat[row0 * N + idx1(k)] -= inv_L;
            }
        }

        double base = p.left_sum[r] + static_cast<double>(p.right_cnt[r]) * r;
        if (p.mid[r]) base += r;
        rhs[row0] = base * inv_L;

        mat[row1 * N + row1] = 1.0;
        if (p.mid[r]) mat[row1 * N + idxAvg1] -= inv_L;
        if (p.right_cnt[r] > 0) {
            for (int k = p.right_k1[r]; k <= p.right_k2[r]; ++k) {
                mat[row1 * N + idx0(k)] -= inv_L;
            }
        }
        rhs[row1] = rhs[row0];
    }

    mat[idxAvg0 * N + idxAvg0] = 1.0;
    for (int r = a; r <= b; ++r) {
        mat[idxAvg0 * N + idx0(r)] -= inv_L;
    }
    rhs[idxAvg0] = 0.0;

    mat[idxAvg1 * N + idxAvg1] = 1.0;
    for (int r = a; r <= b; ++r) {
        mat[idxAvg1 * N + idx1(r)] -= inv_L;
    }
    rhs[idxAvg1] = 0.0;

    vector<double> sol;
    EResult res;
    res.converged = solve_linear(mat, rhs, N, sol);
    res.iterations = 1;
    if (!res.converged) {
        res.value = numeric_limits<double>::quiet_NaN();
        return res;
    }

    double avgA0 = sol[idxAvg0];
    vector<double> A1(M + 1, 0.0);
    for (int r = 1; r <= M; ++r) A1[r] = sol[idx1(r)];

    double total = 0.0;
    for (int ta = n; ta <= m; ++ta) {
        for (int tb = n; tb <= m; ++tb) {
            if (ta < tb) {
                total += ta + A1[tb - ta];
            } else if (tb < ta) {
                total += tb + A1[ta - tb];
            } else {
                total += ta + avgA0;
            }
        }
    }

    res.value = total * inv_L * inv_L;
    return res;
}

double round2(double v) {
    return floor(v * 100.0 + 0.5) / 100.0;
}

bool validate() {
    const double eps = 0.02;
    auto e = compute_E(30, 60);
    if (!isfinite(e.value)) {
        cerr << "Validation failed for E(60,30): non-finite result\n";
        return false;
    }
    double e_round = round2(e.value);
    if (fabs(e_round - 1036.15) > eps) {
        cerr << "Validation failed for E(60,30): got " << fixed << setprecision(2)
             << e_round << '\n';
        return false;
    }

    double s5 = 0.0;
    for (int m = 2; m <= 5; ++m) {
        for (int n = 1; n < m; ++n) {
            auto res = compute_E(n, m);
            if (!isfinite(res.value)) {
                cerr << "Validation failed for S(5): non-finite result\n";
                return false;
            }
            s5 += res.value;
        }
    }
    double s5_round = round2(s5);
    if (fabs(s5_round - 7722.82) > eps) {
        cerr << "Validation failed for S(5): got " << fixed << setprecision(2)
             << s5_round << '\n';
        return false;
    }
    return true;
}

double compute_S(int k, unsigned threads) {
    struct Task {
        int n;
        int m;
    };
    vector<Task> tasks;
    tasks.reserve(k * (k - 1) / 2);
    for (int m = 2; m <= k; ++m) {
        for (int n = 1; n < m; ++n) tasks.push_back({n, m});
    }

    if (threads == 0) threads = 1;
    threads = min<unsigned>(threads, static_cast<unsigned>(tasks.size()));

    atomic<size_t> next{0};
    vector<double> partial(threads, 0.0);
    atomic<bool> all_converged{true};
    atomic<bool> failed{false};
    const size_t chunk = 4;

    auto worker = [&](unsigned idx) {
        double local = 0.0;
        while (true) {
            if (failed.load(memory_order_relaxed)) break;
            size_t start = next.fetch_add(chunk, memory_order_relaxed);
            if (start >= tasks.size()) break;
            size_t end = min(tasks.size(), start + chunk);
            for (size_t i = start; i < end; ++i) {
                const auto& t = tasks[i];
                auto res = compute_E(t.n, t.m);
                if (!isfinite(res.value)) {
                    failed.store(true, memory_order_relaxed);
                    cerr << "Non-finite result for n=" << t.n << " m=" << t.m << '\n';
                    break;
                }
                local += res.value;
                if (!res.converged) all_converged.store(false, memory_order_relaxed);
            }
        }
        partial[idx] = local;
    };

    vector<thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    double total = 0.0;
    for (double v : partial) total += v;

    if (failed.load()) {
        cerr << "Aborting: non-finite intermediate value encountered." << '\n';
        return numeric_limits<double>::quiet_NaN();
    }
    if (!all_converged.load()) {
        cerr << "Warning: solver failed for at least one pair." << '\n';
    }
    return total;
}

} // namespace

int main(int argc, char** argv) {
    if (!validate()) return 1;

    unsigned threads = thread::hardware_concurrency();
    if (argc > 1) {
        int val = atoi(argv[1]);
        if (val > 0) threads = static_cast<unsigned>(val);
    }

    double s100 = compute_S(100, threads);
    cout << fixed << setprecision(2) << s100 << '\n';
    return 0;
}
