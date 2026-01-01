#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace std;

struct Fenwick {
    vector<int> tree;
    int n = 0;

    void reset(int size) {
        n = size;
        tree.assign(n + 1, 0);
    }

    inline void add(int idx, int delta) {
        for (int i = idx + 1; i <= n; i += i & -i) {
            tree[i] += delta;
        }
    }

    inline int sum_prefix(int idx) const {
        int res = 0;
        for (int i = idx + 1; i > 0; i -= i & -i) {
            res += tree[i];
        }
        return res;
    }
};

static int thread_count() {
    unsigned int hc = thread::hardware_concurrency();
    int t = hc == 0 ? 4 : static_cast<int>(hc);
    return min(t, 8);
}

static void build_pos_parallel(long long N, long long step, vector<int> &pos, int threads) {
    if (threads <= 1 || N < 1'000'000) {
        long long r = 0;
        for (int i = 0; i < N; ++i) {
            pos[static_cast<int>(r)] = i;
            r += step;
            if (r >= N) r -= N;
        }
        return;
    }

    vector<thread> pool;
    long long chunk = (N + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
        long long start = t * chunk;
        long long end = min(N, start + chunk);
        pool.emplace_back([&, start, end]() {
            for (long long i = start; i < end; ++i) {
                long long r = (step * i) % N;
                pos[static_cast<int>(r)] = static_cast<int>(i);
            }
        });
    }
    for (auto &th : pool) th.join();
}

static long long count_rect(long long N, int k) {
    long long pow2 = 1LL << k;
    long long T = 2 * N - pow2;
    if (T < 0) return N * N;

    if (T <= N - 1) {
        long long a = N - T - 1;
        long long sum_range = (T * (T + 1)) / 2 + (T + 1) * a;
        long long count_high = N - 1 - T;
        long long sum_high = count_high > 0 ? count_high * N : 0;
        return sum_range + sum_high;
    }

    long long r0 = T - (N - 1);
    if (r0 >= N - 1) return 0;
    long long start = r0 + 1;
    long long end = N - 1;
    if (start > end) return 0;
    long long a = T - N + 1;
    long long count = end - start + 1;
    long long sum_r = (start + end) * count / 2;
    return sum_r - a * count;
}

static long long count_tri_subset(long long N, int k, const vector<int> &pos, vector<int> &d_of_idx) {
    int K = 1 << k;
    int Ksum = K - 2;

    for (int d = 0; d < K; ++d) {
        int idx = pos[static_cast<int>(N - 1 - d)];
        d_of_idx[idx] = d;
    }

    Fenwick bit;
    bit.reset(K);

    long long total = 0;
    int idx_ptr = 0;
    for (int idx = static_cast<int>(N); idx-- > 0;) {
        int d1 = d_of_idx[idx];
        if (d1 < 0) continue;
        int limit = static_cast<int>(N - 1 - idx);
        while (idx_ptr <= limit) {
            int d2 = d_of_idx[idx_ptr];
            if (d2 >= 0) {
                bit.add(d2, 1);
            }
            ++idx_ptr;
        }
        int max_d2 = Ksum - d1 - 1;
        if (max_d2 >= 0) {
            if (max_d2 >= K) max_d2 = K - 1;
            total += bit.sum_prefix(max_d2);
        }
    }

    for (int d = 0; d < K; ++d) {
        int idx = pos[static_cast<int>(N - 1 - d)];
        d_of_idx[idx] = -1;
    }

    return total;
}

static long long count_tri_general(long long N, int k, const vector<int> &pos) {
    long long pow2 = 1LL << k;
    long long T = 2 * N - pow2;
    if (T < 0) {
        return N * (N + 1) / 2;
    }

    Fenwick bit;
    bit.reset(static_cast<int>(N));

    long long total = 0;
    for (int t = static_cast<int>(N - 1); t >= 0; --t) {
        int r_add = t + 1;
        if (r_add < N) {
            bit.add(pos[r_add], 1);
        }
        long long r_q = T - t;
        if (0 <= r_q && r_q < N) {
            int j = pos[static_cast<int>(r_q)];
            int L = static_cast<int>(N - j);
            total += bit.sum_prefix(L - 1);
        }
    }

    if (T < N - 1) {
        long long extra = 0;
        for (int r = static_cast<int>(T + 1); r < N; ++r) {
            extra += static_cast<long long>(N - pos[r]);
        }
        total += extra;
    }

    return total;
}

static long long compute_S_odd(long long N, int threads) {
    if (N % 2 == 0) {
        cerr << "This implementation assumes N is odd.\n";
        exit(1);
    }

    long long total_tri = N * (3 * N + 1) / 2;
    vector<int> pos(static_cast<size_t>(N));
    vector<int> d_of_idx(static_cast<size_t>(N), -1);

    long long S = 0;
    long long F_prev = 0;

    for (int k = 1;; ++k) {
        long long pow2 = 1LL << k;
        if (pow2 > 2 * N) {
            S += total_tri - F_prev;
            break;
        }

        long long step = pow2 % N;
        build_pos_parallel(N, step, pos, threads);

        long long Frect = count_rect(N, k);
        long long Ftri = 0;
        if (pow2 <= N) {
            Ftri = count_tri_subset(N, k, pos, d_of_idx);
        } else {
            Ftri = count_tri_general(N, k, pos);
        }
        long long Fk = Frect + Ftri;

        S += total_tri - F_prev;
        F_prev = Fk;
    }

    return S;
}

static void run_validations() {
    struct Test {
        long long N;
        long long expected;
    };

    vector<Test> tests = {
        {3, 42},
        {5, 126},
        {123, 167178},
        {12345, 3185041956LL},
    };

    for (const auto &t : tests) {
        long long got = compute_S_odd(t.N, 1);
        if (got != t.expected) {
            cerr << "Validation failed for N=" << t.N << ": got " << got
                 << ", expected " << t.expected << "\n";
            exit(1);
        }
    }
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    bool validate = false;
    long long N = 123456789LL;
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        if (s == "--validate") {
            validate = true;
        } else {
            N = stoll(s);
        }
    }

    if (validate) {
        run_validations();
    }

    int threads = thread_count();
    long long ans = compute_S_odd(N, threads);
    cout << ans << "\n";
    return 0;
}
