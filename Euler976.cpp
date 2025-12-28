#include <bits/stdc++.h>
using namespace std;

static constexpr int MOD = 1234567891;

static inline int addmod(int a, int b) {
    long long s = (long long)a + b;
    if (s >= MOD) s -= MOD;
    return (int)s;
}
static inline int submod(int a, int b) {
    long long s = (long long)a - b;
    if (s < 0) s += MOD;
    return (int)s;
}
static inline int mulmod(long long a, long long b) {
    return (int)((a * b) % MOD);
}

static long long modpow(long long a, long long e) {
    long long r = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return r;
}

struct Comb {
    int maxn;
    vector<int> fact;
    vector<int> invfact;

    explicit Comb(int maxn_) : maxn(maxn_), fact(maxn_ + 1), invfact(maxn_ + 1) {
        fact[0] = 1;
        for (int i = 1; i <= maxn; i++) {
            fact[i] = mulmod(fact[i - 1], i);
        }
        invfact[maxn] = (int)modpow(fact[maxn], MOD - 2);
        for (int i = maxn; i >= 1; i--) {
            invfact[i - 1] = mulmod(invfact[i], i);
        }

        if (mulmod(fact[maxn], invfact[maxn]) != 1) {
            cerr << "[FATAL] factorial inversion check failed.\n";
            exit(1);
        }

        for (int n : {0, 1, 2, 10, 123, maxn}) {
            if (n < 0 || n > maxn) continue;
            int c0 = nCk(n, 0);
            int cn = nCk(n, n);
            if (c0 != 1 || cn != 1) {
                cerr << "[FATAL] nCk identity failed at n=" << n << "\n";
                exit(1);
            }
        }

        int inv2 = (MOD + 1) / 2;
        if (mulmod(inv2, 2) != 1) {
            cerr << "[FATAL] inv2 check failed.\n";
            exit(1);
        }
    }

    inline int nCk(int n, int k) const {
        if (k < 0 || k > n) return 0;
        return mulmod(fact[n], mulmod(invfact[k], invfact[n - k]));
    }
};

template <class F>
static int parallel_sum(int L, int R, int num_threads, const F& f) {
    int len = R - L;
    if (len <= 0) return 0;

    if (num_threads <= 1 || len < 100000) {
        long long s = 0;
        for (int i = L; i < R; i++) {
            s += f(i);
            if (s >= (1LL << 62)) s %= MOD;
        }
        return (int)(s % MOD);
    }

    num_threads = min(num_threads, len);
    vector<thread> th;
    vector<int> partial(num_threads, 0);

    for (int t = 0; t < num_threads; t++) {
        int a = L + (long long)len * t / num_threads;
        int b = L + (long long)len * (t + 1) / num_threads;
        th.emplace_back([&, t, a, b]() {
            long long s = 0;
            for (int i = a; i < b; i++) {
                s += f(i);
                if (s >= (1LL << 62)) s %= MOD;
            }
            partial[t] = (int)(s % MOD);
        });
    }
    for (auto& x : th) x.join();

    int ans = 0;
    for (int t = 0; t < num_threads; t++) ans = addmod(ans, partial[t]);
    return ans;
}

static int solve_N_eq_K_div4(int N, int threads) {
    if (N % 4 != 0) {
        cerr << "[FATAL] This solver expects N % 4 == 0.\n";
        exit(1);
    }

    const int No = N / 2;
    const int Na = N / 4;
    const int maxFac = 2 * N;

    Comb C(maxFac);

    const int inv2 = (MOD + 1) / 2;

    int total_even = parallel_sum(0, No + 1, threads, [&](int m) -> int {
        int k = 2 * m;
        int n = N + k - 1;
        return C.nCk(n, k);
    });

    int L1_raw = parallel_sum(0, Na + 1, threads, [&](int r) -> int {
        int a = C.nCk(No, 2 * r);
        int b = C.nCk(3 * No - r, No - r);
        return mulmod(a, b);
    });
    int L1 = mulmod(inv2, L1_raw);

    int L2 = mulmod(inv2, C.nCk(5 * Na, No));

    int lost_even = addmod(L1, L2);

    int Wodd_raw = parallel_sum(0, Na, threads, [&](int r) -> int {
        int a = C.nCk(No, 2 * r + 1);
        int b = C.nCk(3 * No - 1 - r, No - 1 - r);
        return mulmod(a, b);
    });
    int win_odd = mulmod(inv2, Wodd_raw);

    int ans = addmod(submod(total_even, lost_even), win_odd);
    return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int N = 10'000'000;

    int threads = (int)thread::hardware_concurrency();
    if (threads <= 0) threads = 1;

    int ans = solve_N_eq_K_div4(N, threads);
    cout << ans << "\n";
    return 0;
}