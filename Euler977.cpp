#include <bits/stdc++.h>
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

// Precomputed powers for bases 2..N+1, packed into one array.
static vector<int> pow_all;
static vector<int> pow_offset;
static vector<int> pow_len;

static void build_powers(int n) {
    int max_base = n + 1;
    pow_offset.assign(max_base + 1, 0);
    pow_len.assign(max_base + 1, 0);

    long long total_len = 0;
    for (int b = 2; b <= max_base; b++) {
        int max_exp = (n - 1) / (b - 1) + 1;
        pow_len[b] = max_exp + 1;
        pow_offset[b] = (int)total_len;
        total_len += pow_len[b];
    }

    pow_all.assign((size_t)total_len, 0);
    for (int b = 2; b <= max_base; b++) {
        int off = pow_offset[b];
        int len = pow_len[b];
        pow_all[off] = 1;
        for (int e = 1; e < len; e++) {
            pow_all[off + e] = mulmod(pow_all[off + e - 1], b);
        }
    }
}

static inline int pow_base(int base, int exp) {
    if (exp == 0 || base == 1) return 1;
    return pow_all[pow_offset[base] + exp];
}

// Slow O(n^2) solver for validation.
static int solve_slow(int n) {
    long long total = 0;
    for (int t = 0; t < n; t++) {
        int N = n - t;
        for (int l = 1; l <= N; l++) {
            int q = N / l;
            int r = N % l;
            long long P = (modpow(q, l - r) * modpow(q + 1, r)) % MOD;
            if (t == 0) {
                total += P;
            } else {
                total += ((r == 0 ? q - 1 : q) * P) % MOD;
            }
            if (total >= (1LL << 62)) total %= MOD;
        }
    }
    return (int)(total % MOD);
}

// Brute force by enumerating sequences a_1..a_n and checking a_{k+1} = a_{a_k}.
static long long brute_count(int n) {
    vector<int> a(n + 1, 0);
    long long cnt = 0;
    function<void(int)> dfs = [&](int idx) {
        if (idx > n) {
            for (int k = 1; k < n; k++) {
                if (a[k + 1] != a[a[k]]) return;
            }
            cnt++;
            return;
        }
        for (int v = 1; v <= n; v++) {
            a[idx] = v;
            dfs(idx + 1);
        }
    };
    dfs(1);
    return cnt;
}

static int compute_B_range(int n, int l_start, int l_end) {
    long long sum = 0;
    int n1 = n - 1;
    for (int l = l_start; l < l_end; l++) {
        int max_q = n1 / l;
        for (int q = 1; q <= max_q; q++) {
            int N0 = q * l;
            int N1 = (q + 1) * l - 1;
            if (N1 > n1) N1 = n1;
            int R = N1 - N0;

            int powqL = pow_base(q, l);
            int term0 = mulmod(q - 1, powqL);
            int sum_block = term0;

            if (R >= 1) {
                int powqL1 = pow_base(q, l + 1);
                int powqL1minusR = pow_base(q, l + 1 - R);
                int powq1R1 = pow_base(q + 1, R + 1);
                int sum_r = submod(mulmod(powqL1minusR, powq1R1),
                                   mulmod(powqL1, q + 1));
                sum_block = addmod(sum_block, sum_r);
            }

            sum += sum_block;
            if (sum >= (1LL << 62)) sum %= MOD;
        }
    }
    return (int)(sum % MOD);
}

static int solve_fast(int n, int threads) {
    if (n <= 0) return 0;

    int sumA = 0;
    for (int l = 1; l <= n; l++) {
        int q = n / l;
        int r = n % l;
        int P = mulmod(pow_base(q, l - r), pow_base(q + 1, r));
        sumA = addmod(sumA, P);
    }

    if (n == 1) return sumA;

    int n1 = n - 1;
    if (threads <= 1 || n1 < 50'000) {
        int sumB = compute_B_range(n, 1, n);
        return addmod(sumA, sumB);
    }

    threads = min(threads, n1);
    vector<int> partial(threads, 0);
    vector<thread> th;
    th.reserve(threads);

    for (int t = 0; t < threads; t++) {
        int L = 1 + (long long)n1 * t / threads;
        int R = 1 + (long long)n1 * (t + 1) / threads;
        th.emplace_back([&, t, L, R]() {
            partial[t] = compute_B_range(n, L, R);
        });
    }
    for (auto& tt : th) tt.join();

    int sumB = 0;
    for (int v : partial) sumB = addmod(sumB, v);
    return addmod(sumA, sumB);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int N = 1'000'000;

    if (brute_count(7) != 174) {
        cerr << "[FATAL] brute_count(7) validation failed.\n";
        return 1;
    }
    if (solve_slow(100) != 305741269) {
        cerr << "[FATAL] solve_slow(100) validation failed.\n";
        return 1;
    }

    build_powers(N);

    if (solve_fast(7, 1) != 174) {
        cerr << "[FATAL] solve_fast(7) validation failed.\n";
        return 1;
    }

    int threads = (int)thread::hardware_concurrency();
    if (threads <= 0) threads = 1;

    int ans = solve_fast(N, threads);
    cout << ans << "\n";
    return 0;
}
