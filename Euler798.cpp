#include <bits/stdc++.h>
#include <pthread.h>
#include <unistd.h>

using namespace std;

static constexpr int MOD = 1'000'000'007;

static inline int addmod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

static inline int submod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

static long long mod_pow(long long a, long long e) {
    long long r = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return r;
}

static int next_pow2(int x) {
    int p = 1;
    while (p < x) p <<= 1;
    return p;
}

static int get_thread_count() {
    const char* env = getenv("EULER_THREADS");
    if (env && *env) {
        long v = strtol(env, nullptr, 10);
        if (v > 0) return (int)v;
    }
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    return (int)n;
}

static vector<int> compute_distribution_one_suit(int n, int L, bool do_checks) {
    vector<int> arr(L, 0);
    if (n <= 0) return arr;

    if (n == 1) {
        arr[0] = 2;
        return arr;
    }

    int maxInv = n;
    vector<int> inv(maxInv + 1);
    inv[1] = 1;
    for (int i = 2; i <= maxInv; i++) {
        inv[i] = int(MOD - (long long)(MOD / i) * inv[MOD % i] % MOD);
    }
    const int inv2 = (MOD + 1) / 2;

    long long pow2N = mod_pow(2, n - 2);
    arr[0] = int((pow2N + 2) % MOD);

    long long sumDist = 0;
    if (do_checks) sumDist = arr[0];

    long long B = 1;
    long long T = 1;

    int m_max = (n - 1) / 2;

    for (int m = 0; m <= m_max; m++) {
        long long denomA = (long long)n - 1 - 2LL * m;
        long long A = 0;
        if (m >= 1) {
            if (denomA == 0) {
                A = 1;
            } else {
                A = (B * m) % MOD * inv[(int)denomA] % MOD;
            }
        }

        long long P = (B + A) % MOD;
        long long Q = P * ((long long)n - 2LL * m - 1) % MOD * inv[m + 1] % MOD;

        int godd = 2 * m + 1;
        if (godd < n) {
            long long val = (pow2N + Q - T) % MOD;
            if (val < 0) val += MOD;
            arr[godd] = (int)val;
            if (do_checks) sumDist += arr[godd];
        }

        int geven = 2 * m;
        if (m >= 1 && geven < n) {
            long long val = (pow2N + P + B - T) % MOD;
            if (val < 0) val += MOD;
            arr[geven] = (int)val;
            if (do_checks) sumDist += arr[geven];
        }

        if (m == m_max) break;

        long long num1 = (long long)n - 2 - 2LL * m;
        long long num2 = (long long)n - 3 - 2LL * m;
        long long den1 = (long long)m + 1;
        long long den2 = (long long)n - 2 - m;
        B = B * (num1 % MOD) % MOD;
        B = B * (num2 % MOD) % MOD;
        B = B * inv[(int)den1] % MOD;
        B = B * inv[(int)den2] % MOD;

        pow2N = pow2N * inv2 % MOD;

        int m1 = m + 1;
        long long denomA_next = (long long)n - 1 - 2LL * m1;
        long long A_next = 0;
        if (denomA_next == 0) A_next = 1;
        else A_next = (B * m1) % MOD * inv[(int)denomA_next] % MOD;

        T = ((T + A_next) % MOD) * inv2 % MOD;
        T = (T + B) % MOD;
    }

    if (do_checks) {
        sumDist %= MOD;
        long long want = mod_pow(2, n);
        if (sumDist != want) {
            cerr << "[CHECK FAILED] sum(a[g]) != 2^n (mod MOD)\n";
            cerr << "sumDist=" << sumDist << " want=" << want << "\n";
            exit(1);
        }
        if (arr[n - 1] != 1) {
            cerr << "[CHECK FAILED] a[n-1] != 1\n";
            cerr << "a[n-1]=" << arr[n - 1] << "\n";
            exit(1);
        }
    }

    return arr;
}

struct WhtTask {
    vector<int>* a;
    int len;
    int start_block;
    int end_block;
};

static void* wht_worker(void* arg) {
    WhtTask* t = (WhtTask*)arg;
    vector<int>& a = *t->a;
    int len = t->len;
    int step = len << 1;
    for (int b = t->start_block; b < t->end_block; b++) {
        int i = b * step;
        for (int j = 0; j < len; j++) {
            int u = a[i + j];
            int v = a[i + j + len];
            a[i + j] = addmod(u, v);
            a[i + j + len] = submod(u, v);
        }
    }
    return nullptr;
}

static void walsh_hadamard_xor(vector<int>& a) {
    const int n = (int)a.size();
    const int threads = get_thread_count();
    for (int len = 1; len < n; len <<= 1) {
        const int step = len << 1;
        const int blocks = n / step;
        if (threads > 1 && blocks >= 2048) {
            int tcnt = min(threads, blocks);
            vector<pthread_t> th(tcnt);
            vector<WhtTask> tasks(tcnt);
            int per = (blocks + tcnt - 1) / tcnt;
            for (int i = 0; i < tcnt; i++) {
                int start = i * per;
                int end = min(blocks, start + per);
                tasks[i] = WhtTask{&a, len, start, end};
                pthread_create(&th[i], nullptr, wht_worker, &tasks[i]);
            }
            for (int i = 0; i < tcnt; i++) {
                pthread_join(th[i], nullptr);
            }
        } else {
            for (int i = 0; i < n; i += step) {
                for (int j = 0; j < len; j++) {
                    int u = a[i + j];
                    int v = a[i + j + len];
                    a[i + j] = addmod(u, v);
                    a[i + j + len] = submod(u, v);
                }
            }
        }
    }
}

struct PowTask {
    const vector<int>* a;
    int s;
    int start;
    int end;
    long long partial;
};

static void* pow_worker(void* arg) {
    PowTask* t = (PowTask*)arg;
    long long sum = 0;
    for (int i = t->start; i < t->end; i++) {
        sum += mod_pow((*t->a)[i], t->s);
    }
    t->partial = sum % MOD;
    return nullptr;
}

static int solve_C(int n, int s, bool do_checks) {
    int L = next_pow2(n);
    vector<int> A = compute_distribution_one_suit(n, L, do_checks);

    walsh_hadamard_xor(A);

    long long total = 0;
    int threads = get_thread_count();
    if (threads > 1 && L >= 4096) {
        int tcnt = min(threads, L);
        vector<pthread_t> th(tcnt);
        vector<PowTask> tasks(tcnt);
        int per = (L + tcnt - 1) / tcnt;
        for (int i = 0; i < tcnt; i++) {
            int start = i * per;
            int end = min(L, start + per);
            tasks[i] = PowTask{&A, s, start, end, 0};
            pthread_create(&th[i], nullptr, pow_worker, &tasks[i]);
        }
        for (int i = 0; i < tcnt; i++) {
            pthread_join(th[i], nullptr);
            total += tasks[i].partial;
        }
        total %= MOD;
    } else {
        for (int i = 0; i < L; i++) {
            total += mod_pow(A[i], s);
        }
        total %= MOD;
    }

    long long invL = mod_pow(L, MOD - 2);
    return int(total * invL % MOD);
}

static void run_selftest() {
    int a = solve_C(3, 2, true);
    int b = solve_C(13, 4, true);
    if (a != 26) {
        cerr << "[SELFTEST FAILED] C(3,2) expected 26 got " << a << "\n";
        exit(1);
    }
    if (b != 540318329) {
        cerr << "[SELFTEST FAILED] C(13,4) expected 540318329 got " << b << "\n";
        exit(1);
    }
    cerr << "[SELFTEST OK]\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 10'000'000;
    int s = 10'000'000;
    bool do_checks = true;

    if (argc >= 2) {
        string arg1 = argv[1];
        if (arg1 == "--selftest") {
            run_selftest();
            return 0;
        }
    }
    if (argc >= 3) {
        n = stoi(argv[1]);
        s = stoi(argv[2]);
    }

    int ans = solve_C(n, s, do_checks);
    cout << ans << "\n";
    return 0;
}
