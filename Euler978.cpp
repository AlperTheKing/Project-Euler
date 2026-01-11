#include <bits/stdc++.h>
using namespace std;

struct Moments {
    long long second;
    long long third;
};

static Moments compute_moments(int n) {
    long long a0 = 0, a1 = 1;
    long long c0 = 0, c1 = 1;

    if (n == 0) return {a0, c0};
    if (n == 1) return {a1, c1};

    for (int t = 2; t <= n; t++) {
        long long a2 = a1 + a0;
        long long c2 = c1 + 3 * c0;
        a0 = a1;
        a1 = a2;
        c0 = c1;
        c1 = c2;
    }
    return {a1, c1};
}

static long double skewness(int n) {
    Moments m = compute_moments(n);
    long double a = (long double)m.second;
    long double c = (long double)m.third;
    long double var = a - 1.0L;
    if (var <= 0.0L) return 0.0L;
    long double m3 = c - 3.0L * a + 2.0L;
    long double sigma = sqrtl(var);
    return m3 / (sigma * sigma * sigma);
}

static void validate() {
    Moments m5 = compute_moments(5);
    if (m5.second != 5 || m5.third != 19) {
        cerr << "[FATAL] moment recurrence check failed for t=5.\n";
        exit(1);
    }
    long double s5 = skewness(5);
    if (fabsl(s5 - 0.75L) > 1e-12L) {
        cerr << "[FATAL] skewness check failed for t=5.\n";
        exit(1);
    }
    long double s10 = skewness(10);
    if (fabsl(s10 - 2.50997097L) > 5e-9L) {
        cerr << "[FATAL] skewness check failed for t=10.\n";
        exit(1);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    validate();

    const int N = 50;
    long double ans = skewness(N);
    cout << fixed << setprecision(8) << ans << "\n";
    return 0;
}
