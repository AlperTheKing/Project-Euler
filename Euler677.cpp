#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr int MOD = 1'000'000'007;

int mod_pow(long long a, long long e) {
    long long r = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1LL) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1LL;
    }
    return static_cast<int>(r);
}

int mod_norm(long long x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return static_cast<int>(x);
}

void update_square(const vector<int>& s, vector<int>& s_sq, int n, int limit_n) {
    const int s_n = s[n];
    if (s_n == 0) return;
    const int max_k = min(n - 1, limit_n - n);
    for (int k = 1; k <= max_k; ++k) {
        const int v = s[k];
        if (v == 0) continue;
        const int idx = n + k;
        const int add = static_cast<int>((2LL * s_n * v) % MOD);
        int val = s_sq[idx] + add;
        if (val >= MOD) val -= MOD;
        s_sq[idx] = val;
    }
    if (2 * n <= limit_n) {
        const int add = static_cast<int>((1LL * s_n * s_n) % MOD);
        int val = s_sq[2 * n] + add;
        if (val >= MOD) val -= MOD;
        s_sq[2 * n] = val;
    }
}

int conv_cube_at(const vector<int>& s, const vector<int>& s_sq, int t) {
    if (t <= 1) return 0;
    long long sum = 0;
    for (int k = 1; k <= t - 1; ++k) {
        const int a = s[k];
        const int b = s_sq[t - k];
        if (a == 0 || b == 0) continue;
        sum += (1LL * a * b) % MOD;
    }
    return static_cast<int>(sum % MOD);
}

int conv_sc2_at(const vector<int>& s, int t) {
    if (t <= 1) return 0;
    long long sum = 0;
    for (int u = 2; u <= t - 1; u += 2) {
        const int k = t - u;
        if (k < 1) break;
        const int a = s[k];
        const int b = s[u / 2];
        if (a == 0 || b == 0) continue;
        sum += (1LL * a * b) % MOD;
    }
    return static_cast<int>(sum % MOD);
}

vector<int> build_shift(const vector<int>& s, int k, int nmax) {
    vector<int> out(nmax + 1, 0);
    for (int i = 1; i * k <= nmax; ++i) {
        out[i * k] = s[i];
    }
    return out;
}

vector<int> convolution_single(const vector<int>& a, const vector<int>& b, int nmax) {
    vector<int> out(nmax + 1, 0);
    for (int i = 1; i <= nmax; ++i) {
        if (a[i] == 0) continue;
        const int ai = a[i];
        for (int j = 1; j + i <= nmax; ++j) {
            if (b[j] == 0) continue;
            const int idx = i + j;
            const int add = static_cast<int>((1LL * ai * b[j]) % MOD);
            int val = out[idx] + add;
            if (val >= MOD) val -= MOD;
            out[idx] = val;
        }
    }
    return out;
}

vector<int> convolution_parallel(const vector<int>& a, const vector<int>& b, int nmax, unsigned threads) {
    if (threads <= 1 || nmax < 512) {
        return convolution_single(a, b, nmax);
    }
    vector<vector<int>> partials(threads, vector<int>(nmax + 1, 0));
    vector<thread> workers;
    const int chunk = (nmax + static_cast<int>(threads) - 1) / static_cast<int>(threads);
    for (unsigned tid = 0; tid < threads; ++tid) {
        int start = static_cast<int>(tid) * chunk + 1;
        int end = min(nmax, start + chunk - 1);
        workers.emplace_back([&, tid, start, end]() {
            auto& local = partials[tid];
            for (int i = start; i <= end; ++i) {
                if (a[i] == 0) continue;
                const int ai = a[i];
                for (int j = 1; j + i <= nmax; ++j) {
                    if (b[j] == 0) continue;
                    const int idx = i + j;
                    const int add = static_cast<int>((1LL * ai * b[j]) % MOD);
                    int val = local[idx] + add;
                    if (val >= MOD) val -= MOD;
                    local[idx] = val;
                }
            }
        });
    }
    for (auto& th : workers) th.join();
    vector<int> out(nmax + 1, 0);
    for (unsigned tid = 0; tid < threads; ++tid) {
        const auto& local = partials[tid];
        for (int i = 1; i <= nmax; ++i) {
            int val = out[i] + local[i];
            if (val >= MOD) val -= MOD;
            out[i] = val;
        }
    }
    return out;
}

vector<int> convolution_stride(const vector<int>& a, const vector<int>& base, int stride, int nmax, unsigned threads) {
    const int max_k = nmax / stride;
    if (threads <= 1 || max_k < 256) {
        vector<int> out(nmax + 1, 0);
        for (int k = 1; k <= max_k; ++k) {
            const int coeff = base[k];
            if (coeff == 0) continue;
            const int offset = stride * k;
            for (int j = 1; j + offset <= nmax; ++j) {
                const int aj = a[j];
                if (aj == 0) continue;
                const int idx = j + offset;
                const int add = static_cast<int>((1LL * aj * coeff) % MOD);
                int val = out[idx] + add;
                if (val >= MOD) val -= MOD;
                out[idx] = val;
            }
        }
        return out;
    }
    vector<vector<int>> partials(threads, vector<int>(nmax + 1, 0));
    vector<thread> workers;
    const int chunk = (max_k + static_cast<int>(threads) - 1) / static_cast<int>(threads);
    for (unsigned tid = 0; tid < threads; ++tid) {
        int start = static_cast<int>(tid) * chunk + 1;
        int end = min(max_k, start + chunk - 1);
        workers.emplace_back([&, tid, start, end]() {
            auto& local = partials[tid];
            for (int k = start; k <= end; ++k) {
                const int coeff = base[k];
                if (coeff == 0) continue;
                const int offset = stride * k;
                for (int j = 1; j + offset <= nmax; ++j) {
                    const int aj = a[j];
                    if (aj == 0) continue;
                    const int idx = j + offset;
                    const int add = static_cast<int>((1LL * aj * coeff) % MOD);
                    int val = local[idx] + add;
                    if (val >= MOD) val -= MOD;
                    local[idx] = val;
                }
            }
        });
    }
    for (auto& th : workers) th.join();
    vector<int> out(nmax + 1, 0);
    for (unsigned tid = 0; tid < threads; ++tid) {
        const auto& local = partials[tid];
        for (int i = 1; i <= nmax; ++i) {
            int val = out[i] + local[i];
            if (val >= MOD) val -= MOD;
            out[i] = val;
        }
    }
    return out;
}

struct PlantedData {
    vector<int> p_r;
    vector<int> p_b;
    vector<int> p_y;
    vector<int> s_all;
    vector<int> s_no_y;
    vector<int> s_all_sq;
    vector<int> s_no_y_sq;
};

PlantedData build_planted(int nmax) {
    const int inv2 = mod_pow(2, MOD - 2);
    const int inv6 = mod_pow(6, MOD - 2);
    vector<int> p_r(nmax + 1, 0);
    vector<int> p_b(nmax + 1, 0);
    vector<int> p_y(nmax + 1, 0);
    vector<int> s_all(nmax + 1, 0);
    vector<int> s_no_y(nmax + 1, 0);
    vector<int> s_all_sq(nmax + 1, 0);
    vector<int> s_no_y_sq(nmax + 1, 0);

    for (int n = 1; n <= nmax; ++n) {
        const int t = n - 1;
        const int base = (t == 0) ? 1 : 0;

        const int s1_all = s_all[t];
        const int s2_all = s_all_sq[t];
        const int s3_all = conv_cube_at(s_all, s_all_sq, t);
        const int sc2_all = conv_sc2_at(s_all, t);
        const int c2_all = (t % 2 == 0) ? s_all[t / 2] : 0;
        const int c3_all = (t % 3 == 0) ? s_all[t / 3] : 0;

        const int m1_all = s1_all;
        const int m2_all = static_cast<int>((1LL * (s2_all + c2_all) % MOD) * inv2 % MOD);
        const int m3_all = static_cast<int>((1LL * (s3_all + 3LL * sc2_all + 2LL * c3_all) % MOD) * inv6 % MOD);

        const int s1_no = s_no_y[t];
        const int s2_no = s_no_y_sq[t];
        const int s3_no = conv_cube_at(s_no_y, s_no_y_sq, t);
        const int sc2_no = conv_sc2_at(s_no_y, t);
        const int c2_no = (t % 2 == 0) ? s_no_y[t / 2] : 0;
        const int c3_no = (t % 3 == 0) ? s_no_y[t / 3] : 0;

        const int m1_no = s1_no;
        const int m2_no = static_cast<int>((1LL * (s2_no + c2_no) % MOD) * inv2 % MOD);
        const int m3_no = static_cast<int>((1LL * (s3_no + 3LL * sc2_no + 2LL * c3_no) % MOD) * inv6 % MOD);

        p_r[n] = mod_norm(static_cast<long long>(base) + m1_all + m2_all + m3_all);
        p_b[n] = mod_norm(static_cast<long long>(base) + m1_all + m2_all);
        p_y[n] = mod_norm(static_cast<long long>(base) + m1_no + m2_no);

        s_all[n] = mod_norm(static_cast<long long>(p_r[n]) + p_b[n] + p_y[n]);
        s_no_y[n] = mod_norm(static_cast<long long>(p_r[n]) + p_b[n]);

        update_square(s_all, s_all_sq, n, nmax);
        update_square(s_no_y, s_no_y_sq, n, nmax);
    }

    PlantedData out;
    out.p_r = std::move(p_r);
    out.p_b = std::move(p_b);
    out.p_y = std::move(p_y);
    out.s_all = std::move(s_all);
    out.s_no_y = std::move(s_no_y);
    out.s_all_sq = std::move(s_all_sq);
    out.s_no_y_sq = std::move(s_no_y_sq);
    return out;
}

}  // namespace

int main() {
    const int nmax = 10000;
    const int inv2 = mod_pow(2, MOD - 2);
    const int inv6 = mod_pow(6, MOD - 2);
    const int inv24 = mod_pow(24, MOD - 2);
    unsigned threads = thread::hardware_concurrency();
    if (threads == 0) threads = 1;

    PlantedData planted = build_planted(nmax);
    const vector<int>& p_all = planted.s_all;
    const vector<int>& p_no_y = planted.s_no_y;
    const vector<int>& p_y = planted.p_y;

    vector<int> c2_all = build_shift(p_all, 2, nmax);
    vector<int> c3_all = build_shift(p_all, 3, nmax);
    vector<int> c4_all = build_shift(p_all, 4, nmax);

    vector<int> c2_no = build_shift(p_no_y, 2, nmax);
    vector<int> c3_no = build_shift(p_no_y, 3, nmax);

    vector<int> all_cube = convolution_parallel(planted.s_all_sq, p_all, nmax, threads);
    vector<int> all_four = convolution_parallel(planted.s_all_sq, planted.s_all_sq, nmax, threads);
    vector<int> all_c2 = convolution_stride(p_all, p_all, 2, nmax, threads);
    vector<int> all_sq_c2 = convolution_stride(planted.s_all_sq, p_all, 2, nmax, threads);
    vector<int> all_c3 = convolution_stride(p_all, p_all, 3, nmax, threads);

    vector<int> no_cube = convolution_parallel(planted.s_no_y_sq, p_no_y, nmax, threads);
    vector<int> no_c2 = convolution_stride(p_no_y, p_no_y, 2, nmax, threads);

    vector<int> c2_sq_all(nmax + 1, 0);
    for (int i = 1; 2 * i <= nmax; ++i) {
        c2_sq_all[2 * i] = planted.s_all_sq[i];
    }

    vector<int> m1_all = p_all;
    vector<int> m2_all(nmax + 1, 0);
    vector<int> m3_all(nmax + 1, 0);
    vector<int> m4_all(nmax + 1, 0);

    vector<int> m1_no = p_no_y;
    vector<int> m2_no(nmax + 1, 0);
    vector<int> m3_no(nmax + 1, 0);

    for (int t = 0; t <= nmax; ++t) {
        m2_all[t] = static_cast<int>((1LL * (planted.s_all_sq[t] + c2_all[t]) % MOD) * inv2 % MOD);
        m2_no[t] = static_cast<int>((1LL * (planted.s_no_y_sq[t] + c2_no[t]) % MOD) * inv2 % MOD);

        m3_all[t] = static_cast<int>((1LL * (all_cube[t] + 3LL * all_c2[t] + 2LL * c3_all[t]) % MOD) * inv6 % MOD);
        m3_no[t] = static_cast<int>((1LL * (no_cube[t] + 3LL * no_c2[t] + 2LL * c3_no[t]) % MOD) * inv6 % MOD);

        m4_all[t] = static_cast<int>((1LL * (all_four[t] +
                                             6LL * all_sq_c2[t] +
                                             3LL * c2_sq_all[t] +
                                             8LL * all_c3[t] +
                                             6LL * c4_all[t]) % MOD) * inv24 % MOD);
    }

    vector<int> a_all(nmax + 1, 0);
    vector<int> a_r(nmax + 1, 0);
    vector<int> a_b(nmax + 1, 0);
    vector<int> a_y(nmax + 1, 0);

    for (int n = 1; n <= nmax; ++n) {
        const int t = n - 1;
        const int base = (t == 0) ? 1 : 0;
        a_r[n] = mod_norm(static_cast<long long>(base) + m1_all[t] + m2_all[t] + m3_all[t] + m4_all[t]);
        a_b[n] = mod_norm(static_cast<long long>(base) + m1_all[t] + m2_all[t] + m3_all[t]);
        a_y[n] = mod_norm(static_cast<long long>(base) + m1_no[t] + m2_no[t] + m3_no[t]);
        a_all[n] = mod_norm(static_cast<long long>(a_r[n]) + a_b[n] + a_y[n]);
    }

    vector<int> p_y_sq = convolution_parallel(p_y, p_y, nmax, threads);
    vector<int> p_all_x2 = build_shift(p_all, 2, nmax);
    vector<int> p_y_x2 = build_shift(p_y, 2, nmax);

    vector<int> g(nmax + 1, 0);
    for (int n = 1; n <= nmax; ++n) {
        const int d = mod_norm(static_cast<long long>(planted.s_all_sq[n]) - p_y_sq[n]);
        const int e = mod_norm(static_cast<long long>(planted.s_all_sq[n]) + p_all_x2[n] - p_y_sq[n] - p_y_x2[n]);
        const int e_half = static_cast<int>((1LL * e) * inv2 % MOD);
        g[n] = mod_norm(static_cast<long long>(a_all[n]) + e_half - d);
    }

    struct Check {
        int n;
        int expected;
    };
    const Check checks[] = {
        {2, 5},
        {3, 15},
        {4, 57},
        {10, 710249},
        {100, 919747298},
    };
    for (const auto& chk : checks) {
        if (g[chk.n] != chk.expected) {
            cerr << "Validation failed for g(" << chk.n << "): got "
                 << g[chk.n] << " expected " << chk.expected << "\n";
            return 1;
        }
    }

    cout << g[nmax] << "\n";
    return 0;
}
