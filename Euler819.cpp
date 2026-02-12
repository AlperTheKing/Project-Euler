#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

static long double expected_steps(int n) {
    std::vector<std::vector<long double>> p(n + 1, std::vector<long double>(n + 1, 0.0L));
    p[0][0] = 1.0L;

    const long double N = static_cast<long double>(n);

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j <= k; ++j) {
            const long double cur = p[k][j];
            if (cur == 0.0L) {
                continue;
            }
            p[k + 1][j] += cur * (static_cast<long double>(j) / N);
            p[k + 1][j + 1] += cur * (static_cast<long double>(n - j) / N);
        }
    }

    std::vector<long double> e(n + 1, 0.0L);
    e[1] = 0.0L;

    for (int k = 2; k <= n; ++k) {
        long double num = 1.0L;
        for (int j = 1; j < k; ++j) {
            num += p[k][j] * e[j];
        }
        const long double den = 1.0L - p[k][k];
        e[k] = num / den;
    }

    return e[n];
}

int main() {
    const long double e3 = expected_steps(3);
    assert(std::fabsl(e3 - (27.0L / 7.0L)) < 1e-14L);

    const long double e5 = expected_steps(5);
    assert(std::fabsl(e5 - (468125.0L / 60701.0L)) < 1e-12L);

    const long double ans = expected_steps(1000);
    std::cout << std::fixed << std::setprecision(6) << ans << '\n';
    return 0;
}
