#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

struct Trans {
    int nxt;
    long double p;
};

static long double R_value(int n) {
    std::vector<Trans> t5(n), t6(n);
    for (int r = 1; r < n; ++r) {
        int m5 = 5 * r;
        int n5 = m5 % n;
        t5[r] = {n5, static_cast<long double>(n5) / static_cast<long double>(m5)};

        int m6 = 6 * r;
        int n6 = m6 % n;
        t6[r] = {n6, static_cast<long double>(n6) / static_cast<long double>(m6)};
    }

    std::vector<long double> v(n, 0.0L), nv(n, 0.0L);

    for (int iter = 0; iter < 2000; ++iter) {
        long double diff = 0.0L;
        for (int r = 1; r < n; ++r) {
            long double c5 = 1.0L + t5[r].p * v[t5[r].nxt];
            long double c6 = 1.0L + t6[r].p * v[t6[r].nxt];
            nv[r] = (c5 < c6) ? c5 : c6;
            long double d = std::fabsl(nv[r] - v[r]);
            if (d > diff) diff = d;
        }
        v.swap(nv);
        if (diff < 1e-18L) break;
    }

    return v[1];
}

static long double S_value(int n) {
    long double s = 0.0L;
    for (int k = 2; k <= n; ++k) s += R_value(k);
    return s;
}

int main() {
    assert(std::fabsl(R_value(8) - 2.083333333333333333L) < 1e-12L);
    assert(std::fabsl(R_value(28) - 2.142476190476190476L) < 1e-12L);
    assert(std::fabsl(S_value(30) - 56.054622093376053963L) < 1e-9L);

    long double ans = S_value(1000);
    std::cout << std::fixed << std::setprecision(6) << ans << '\n';
    return 0;
}
