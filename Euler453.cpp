#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;

using int64 = long long;
using i128 = __int128_t;

static constexpr int64 MOD = 135707531;

static int64 mod_add(int64 a, int64 b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}

static int64 mod_sub(int64 a, int64 b) {
    a -= b;
    if (a < 0) a += MOD;
    return a;
}

static int64 mod_mul(int64 a, int64 b) {
    return static_cast<int64>((static_cast<i128>(a) * b) % MOD);
}

static int64 mod_pow(int64 base, int64 exp) {
    int64 res = 1 % MOD;
    int64 cur = base % MOD;
    while (exp > 0) {
        if (exp & 1) res = mod_mul(res, cur);
        cur = mod_mul(cur, cur);
        exp >>= 1;
    }
    return res;
}

static int64 mod_inv(int64 x) {
    return mod_pow(x, MOD - 2);
}

static string to_string_i128(i128 v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    if (neg) v = -v;
    string s;
    while (v > 0) {
        int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }
    if (neg) s.push_back('-');
    reverse(s.begin(), s.end());
    return s;
}

struct ModStats {
    int64 c3 = 0;
    int64 c4 = 0;
    int64 sumg = 0;
    int64 sum_area2 = 0; // Sum of double areas (mod MOD)
};

struct ExactStats {
    i128 c3 = 0;
    i128 c4 = 0;
    i128 sumg = 0;
    i128 sum_area2 = 0; // Sum of double areas (exact)
};

static ModStats compute_stats_mod(int m, int n, int threads) {
    ModStats total;

    const int max_g = max(m, n);
    vector<int64> g_mod(max_g + 1, 0);
    vector<int64> g1_mod(max_g + 1, 0);
    vector<int64> g2_mod(max_g + 1, 0);
    vector<int64> comb2_mod(max_g + 1, 0);
    for (int g = 0; g <= max_g; g++) {
        g_mod[g] = g % MOD;
        g1_mod[g] = (g >= 1) ? (g - 1) % MOD : 0;
        g2_mod[g] = mod_mul(g_mod[g], g_mod[g]);
        if (g >= 2) {
            int64 comb = static_cast<int64>(g - 1) * static_cast<int64>(g - 2) / 2;
            comb2_mod[g] = comb % MOD;
        }
    }

    vector<int64> mx_mod(m + 1, 0);
    vector<int64> ny_mod(n + 1, 0);
    vector<int64> dx2_mod(m + 1, 0);
    vector<int64> dy2_mod(n + 1, 0);

    for (int dx = 0; dx <= m; dx++) {
        mx_mod[dx] = (static_cast<int64>(m + 1 - dx)) % MOD;
        dx2_mod[dx] = mod_mul(dx % MOD, dx % MOD);
    }
    for (int dy = 0; dy <= n; dy++) {
        ny_mod[dy] = (static_cast<int64>(n + 1 - dy)) % MOD;
        dy2_mod[dy] = mod_mul(dy % MOD, dy % MOD);
    }

    const int64 n_plus_1_mod = (n + 1) % MOD;
    const int64 m_plus_1_mod = (m + 1) % MOD;

    // dx > 0, dy = 0
    for (int dx = 1; dx <= m; dx++) {
        int g = dx;
        int64 pairs = mod_mul(n_plus_1_mod, mx_mod[dx]);
        total.sumg = mod_add(total.sumg, mod_mul(pairs, g_mod[g]));
        if (g >= 2) total.c3 = mod_add(total.c3, mod_mul(pairs, g1_mod[g]));
        if (g >= 3) total.c4 = mod_add(total.c4, mod_mul(pairs, comb2_mod[g]));
    }

    // dx = 0, dy > 0
    for (int dy = 1; dy <= n; dy++) {
        int g = dy;
        int64 pairs = mod_mul(m_plus_1_mod, ny_mod[dy]);
        total.sumg = mod_add(total.sumg, mod_mul(pairs, g_mod[g]));
        if (g >= 2) total.c3 = mod_add(total.c3, mod_mul(pairs, g1_mod[g]));
        if (g >= 3) total.c4 = mod_add(total.c4, mod_mul(pairs, comb2_mod[g]));
    }

    int T = threads;
    if (T <= 0) T = static_cast<int>(thread::hardware_concurrency());
    if (T <= 0) T = 4;
    T = min(T, m);
    vector<ModStats> locals(T);
    vector<thread> pool;

    auto worker = [&](int tid, int start_dx, int end_dx) {
        ModStats local;
        for (int dx = start_dx; dx <= end_dx; dx++) {
            int64 mx = mx_mod[dx];
            int64 dx2 = dx2_mod[dx];
            for (int dy = 1; dy <= n; dy++) {
                int g = std::gcd(dx, dy);
                int64 w = mod_mul(mx, ny_mod[dy]);
                int64 pairs = mod_add(w, w); // 2*w

                local.sumg = mod_add(local.sumg, mod_mul(pairs, g_mod[g]));
                if (g >= 2) local.c3 = mod_add(local.c3, mod_mul(pairs, g1_mod[g]));
                if (g >= 3) local.c4 = mod_add(local.c4, mod_mul(pairs, comb2_mod[g]));

                // Sum of double areas uses: (1/3) * sum w*(i^2 + j^2 + 11*i^2*j^2 - gcd^2).
                int64 dy2 = dy2_mod[dy];
                int64 term = mod_add(dx2, dy2);
                term = mod_add(term, mod_mul(11, mod_mul(dx2, dy2)));
                term = mod_sub(term, g2_mod[g]);
                local.sum_area2 = mod_add(local.sum_area2, mod_mul(w, term));
            }
        }
        locals[tid] = local;
    };

    int chunk = (m + T - 1) / T;
    for (int t = 0; t < T; t++) {
        int start_dx = t * chunk + 1;
        int end_dx = min(m, (t + 1) * chunk);
        if (start_dx > end_dx) continue;
        pool.emplace_back(worker, t, start_dx, end_dx);
    }
    for (auto& th : pool) th.join();

    for (const auto& st : locals) {
        total.c3 = mod_add(total.c3, st.c3);
        total.c4 = mod_add(total.c4, st.c4);
        total.sumg = mod_add(total.sumg, st.sumg);
        total.sum_area2 = mod_add(total.sum_area2, st.sum_area2);
    }

    total.sum_area2 = mod_mul(total.sum_area2, mod_inv(3));
    return total;
}

static ExactStats compute_stats_exact(int m, int n) {
    ExactStats total;

    // dx > 0, dy = 0
    for (int dx = 1; dx <= m; dx++) {
        i128 pairs = static_cast<i128>(n + 1) * (m + 1 - dx);
        int g = dx;
        total.sumg += pairs * g;
        if (g >= 2) total.c3 += pairs * (g - 1);
        if (g >= 3) total.c4 += pairs * (g - 1) * (g - 2) / 2;
    }

    // dx = 0, dy > 0
    for (int dy = 1; dy <= n; dy++) {
        i128 pairs = static_cast<i128>(m + 1) * (n + 1 - dy);
        int g = dy;
        total.sumg += pairs * g;
        if (g >= 2) total.c3 += pairs * (g - 1);
        if (g >= 3) total.c4 += pairs * (g - 1) * (g - 2) / 2;
    }

    for (int dx = 1; dx <= m; dx++) {
        for (int dy = 1; dy <= n; dy++) {
            i128 w = static_cast<i128>(m + 1 - dx) * (n + 1 - dy);
            int g = std::gcd(dx, dy);
            i128 pairs = 2 * w;

            total.sumg += pairs * g;
            if (g >= 2) total.c3 += pairs * (g - 1);
            if (g >= 3) total.c4 += pairs * (g - 1) * (g - 2) / 2;

            i128 dx2 = static_cast<i128>(dx) * dx;
            i128 dy2 = static_cast<i128>(dy) * dy;
            i128 term = dx2 + dy2 + 11 * dx2 * dy2 - static_cast<i128>(g) * g;
            total.sum_area2 += w * term;
        }
    }

    total.sum_area2 /= 3;
    return total;
}

static int64 nC2_mod(int64 n) {
    return mod_mul(mod_mul(n, mod_sub(n, 1)), mod_inv(2));
}

static int64 nC3_mod(int64 n) {
    return mod_mul(mod_mul(mod_mul(n, mod_sub(n, 1)), mod_sub(n, 2)), mod_inv(6));
}

static int64 nC4_mod(int64 n) {
    int64 v = mod_mul(n, mod_sub(n, 1));
    v = mod_mul(v, mod_sub(n, 2));
    v = mod_mul(v, mod_sub(n, 3));
    return mod_mul(v, mod_inv(24));
}

static int64 compute_q_mod(int m, int n, int threads) {
    ModStats stats = compute_stats_mod(m, n, threads);

    int64 N = mod_mul((m + 1) % MOD, (n + 1) % MOD);
    int64 C2 = nC2_mod(N);
    int64 C3n = nC3_mod(N);
    int64 C4n = nC4_mod(N);

    int64 S_coll = mod_sub(mod_mul(mod_sub(N, 3), stats.c3), mod_mul(3, stats.c4));

    // LineSum = sum k*C(k+1,3) simplifies to 2*C2 + 6*C3 + 4*C4.
    int64 line_sum = mod_add(mod_mul(2, C2), mod_add(mod_mul(6, stats.c3), mod_mul(4, stats.c4)));

    int64 sumB = mod_sub(mod_mul(N, stats.sumg), line_sum);

    int64 ans = C4n;
    ans = mod_sub(ans, S_coll);
    ans = mod_add(ans, mod_mul(2, C3n));
    ans = mod_sub(ans, mod_mul(2, stats.c3));
    ans = mod_add(ans, stats.sum_area2);
    ans = mod_sub(ans, sumB);
    return ans;
}

static i128 nC2_exact(i128 n) { return n * (n - 1) / 2; }
static i128 nC3_exact(i128 n) { return n * (n - 1) * (n - 2) / 6; }
static i128 nC4_exact(i128 n) { return n * (n - 1) * (n - 2) * (n - 3) / 24; }

static i128 compute_q_exact(int m, int n) {
    ExactStats stats = compute_stats_exact(m, n);
    i128 N = static_cast<i128>(m + 1) * (n + 1);
    i128 C2 = nC2_exact(N);
    i128 C3n = nC3_exact(N);
    i128 C4n = nC4_exact(N);

    i128 S_coll = (N - 3) * stats.c3 - 3 * stats.c4;
    i128 line_sum = 2 * C2 + 6 * stats.c3 + 4 * stats.c4;
    i128 sumB = N * stats.sumg - line_sum;

    i128 ans = C4n - S_coll + 2 * C3n - 2 * stats.c3 + stats.sum_area2 - sumB;
    return ans;
}

static void run_checks() {
    struct Test {
        int m;
        int n;
        i128 expected;
    };
    vector<Test> tests = {
        {2, 2, 94},
        {3, 7, 39590},
        {12, 3, 309000},
        {123, 45, static_cast<i128>(70542215894646LL)}
    };

    cout << "Running validation checks...\n";
    for (const auto& t : tests) {
        i128 got = compute_q_exact(t.m, t.n);
        cout << "Q(" << t.m << ", " << t.n << ") = " << to_string_i128(got)
             << " (Expected: " << to_string_i128(t.expected) << ")";
        if (got == t.expected) {
            cout << " [PASS]\n";
        } else {
            cout << " [FAIL]\n";
        }
        assert(got == t.expected);
    }
    cout << "All validation checks passed!\n";
    cout << "------------------------------------------------\n";
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    run_checks();

    const int m = 12345;
    const int n = 6789;

    int threads = static_cast<int>(thread::hardware_concurrency());
    int64 result = compute_q_mod(m, n, threads);

    cout << "Q(" << m << ", " << n << ") mod " << MOD << " = " << result << "\n";
    return 0;
}
