#include <array>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Project Euler 603: Substring Sums
// Use the recurrence for substrings ending at j:
//   F_j = 10*F_{j-1} + d_j * j, and S = sum_j F_j.
// Model the update on (F,S,j,1) with a 4x4 matrix; build the matrix for P(10^6) and
// exponentiate it to apply 10^12 repetitions modulo 1e9+7.

using u64 = std::uint64_t;

static constexpr int MOD = 1000000007;

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

struct Mat {
    int a[4][4];
};

static Mat mat_identity() {
    Mat m{};
    for (int i = 0; i < 4; ++i) m.a[i][i] = 1;
    return m;
}

static Mat mat_mul(const Mat& x, const Mat& y) {
    Mat z{};
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            const int xik = x.a[i][k];
            if (!xik) continue;
            for (int j = 0; j < 4; ++j) {
                z.a[i][j] = addmod(z.a[i][j], mulmod(xik, y.a[k][j]));
            }
        }
    }
    return z;
}

static std::array<int, 4> mat_vec_mul(const Mat& m, const std::array<int, 4>& v) {
    std::array<int, 4> r{0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {
        long long acc = 0;
        for (int k = 0; k < 4; ++k) acc += (long long)m.a[i][k] * v[k];
        r[i] = (int)(acc % MOD);
    }
    return r;
}

static void feed_digit(Mat& m, int d) {
    // Left-multiply the current matrix by the per-digit transform.
    // Update columns using the state update rules.
    for (int col = 0; col < 4; ++col) {
        const int f = m.a[0][col];
        const int s = m.a[1][col];
        const int t = m.a[2][col];
        const int one = m.a[3][col];

        const int dt = mulmod(d, t);
        const int done = mulmod(d, one);
        const int tenf = mulmod(10, f);

        const int f2 = addmod(addmod(tenf, dt), done);
        const int s2 = addmod(s, f2);
        const int t2 = addmod(t, one);

        m.a[0][col] = f2;
        m.a[1][col] = s2;
        m.a[2][col] = t2;
        m.a[3][col] = one;
    }
}

static void feed_number(Mat& m, int x) {
    // Feed decimal digits of x in most-significant-first order, no allocations.
    int buf[16];
    int len = 0;
    while (x > 0) {
        buf[len++] = x % 10;
        x /= 10;
    }
    for (int i = len - 1; i >= 0; --i) feed_digit(m, buf[i]);
}

static Mat matrix_from_string(const std::string& s) {
    Mat m = mat_identity();
    for (char ch : s) feed_digit(m, ch - '0');
    return m;
}

static int substring_sum_mod(const std::string& s) {
    long long f = 0;
    long long ans = 0;
    for (int i = 0; i < (int)s.size(); ++i) {
        const int d = s[i] - '0';
        const long long j = i + 1; // 1-indexed
        f = (10 * f + (long long)d * j) % MOD;
        ans += f;
        if (ans >= MOD) ans %= MOD;
    }
    return (int)(ans % MOD);
}

static int upper_bound_nth_prime(int n) {
    // For n >= 6: p_n < n (log n + log log n) (Rosser-Schoenfeld).
    if (n < 6) return 15;
    const double nd = (double)n;
    const double bound = nd * (std::log(nd) + std::log(std::log(nd)));
    return (int)std::ceil(bound + 100.0); // generous margin for rounding
}

static Mat matrix_for_P_primes(int n) {
    // Build the per-block matrix for P(n): concatenation of the first n primes.
    int limit = upper_bound_nth_prime(n);
    for (;;) {
        std::vector<std::uint8_t> comp((limit >> 1) + 1, 0); // only odds; idx = x>>1
        const int r = (int)std::sqrt((double)limit);
        for (int p = 3; p <= r; p += 2) {
            if (comp[p >> 1]) continue;
            const int step = p << 1;
            for (int x = p * p; x <= limit; x += step) comp[x >> 1] = 1;
        }

        Mat m = mat_identity();
        int count = 0;

        // prime 2
        ++count;
        feed_number(m, 2);
        if (count == n) return m;

        for (int x = 3; x <= limit; x += 2) {
            if (comp[x >> 1]) continue;
            ++count;
            feed_number(m, x);
            if (count == n) return m;
        }

        // Not enough primes; increase the limit and retry (shouldn't happen with the bound).
        limit = (int)(limit * 1.2) + 1000;
    }
}

int main() {
    // Validations from the statement / definition.
    if (substring_sum_mod("2024") != 2304) {
        std::cerr << "Validation failed: S(2024)\n";
        return 1;
    }

    // Cross-check matrix method on a small case: P(7) repeated 3 times.
    const std::string p7 = "2357111317";
    const std::string c73 = p7 + p7 + p7;
    const int naive = substring_sum_mod(c73);
    Mat m7 = matrix_from_string(p7);
    std::array<int, 4> v{0, 0, 0, 1};
    u64 k = 3;
    while (k) {
        if (k & 1) v = mat_vec_mul(m7, v);
        m7 = mat_mul(m7, m7);
        k >>= 1;
    }
    if (v[1] != naive) {
        std::cerr << "Validation failed: matrix repetition\n";
        return 1;
    }

    const int N = 1000000;
    const u64 K = 1000000000000ULL;

    Mat base = matrix_for_P_primes(N);
    std::array<int, 4> state{0, 0, 0, 1};

    u64 e = K;
    while (e) {
        if (e & 1) state = mat_vec_mul(base, state);
        base = mat_mul(base, base);
        e >>= 1;
    }

    std::cout << state[1] << "\n";
    return 0;
}
