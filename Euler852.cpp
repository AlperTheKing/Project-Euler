#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using Real = long double;

static Real H_value(Real p0, int depth_limit) {
    if (p0 <= 0.0L || p0 >= 1.0L) return 20.0L;

    const Real ln2 = std::log(2.0L);
    const Real ln3 = std::log(3.0L);
    const Real log_odds0 = std::log(p0 / (1.0L - p0));

    std::vector<Real> next(depth_limit + 1);
    int n = depth_limit;
    Real base = log_odds0 - n * ln2;
    for (int h = 0; h <= n; ++h) {
        Real z = base + h * ln3;
        Real p;
        if (z > 50.0L) p = 1.0L;
        else if (z < -50.0L) p = 0.0L;
        else {
            Real o = std::expl(z);
            p = o / (1.0L + o);
        }
        next[h] = 70.0L * std::max(p, 1.0L - p) - 50.0L;
    }

    for (n = depth_limit - 1; n >= 0; --n) {
        std::vector<Real> cur(n + 1);
        Real base_n = log_odds0 - n * ln2;
        for (int h = 0; h <= n; ++h) {
            Real z = base_n + h * ln3;
            Real p;
            if (z > 50.0L) p = 1.0L;
            else if (z < -50.0L) p = 0.0L;
            else {
                Real o = std::expl(z);
                p = o / (1.0L + o);
            }

            Real stop = 70.0L * std::max(p, 1.0L - p) - 50.0L;
            Real qh = 0.5L + 0.25L * p;
            Real cont = -1.0L + qh * next[h + 1] + (1.0L - qh) * next[h];
            cur[h] = std::max(stop, cont);
        }
        next.swap(cur);
    }

    return next[0];
}

static Real S_value(int N, int depth_limit) {
    std::vector<std::vector<Real>> H(N + 1, std::vector<Real>(N + 1, 20.0L));
    for (int u = 1; u <= N; ++u) {
        for (int f = 1; f <= N; ++f) {
            H[u][f] = H_value(static_cast<Real>(u) / static_cast<Real>(u + f), depth_limit);
        }
    }

    std::vector<std::vector<Real>> V(N + 1, std::vector<Real>(N + 1, 0.0L));
    for (int u = 1; u <= N; ++u) V[u][0] = 20.0L * u;
    for (int f = 1; f <= N; ++f) V[0][f] = 20.0L * f;

    for (int u = 1; u <= N; ++u) {
        for (int f = 1; f <= N; ++f) {
            Real p = static_cast<Real>(u) / static_cast<Real>(u + f);
            V[u][f] = p * V[u - 1][f] + (1.0L - p) * V[u][f - 1] + H[u][f];
        }
    }

    return V[N][N];
}

int main() {
    const int depth_limit = 220;

    Real s1 = S_value(1, depth_limit);
    assert(std::fabsl(s1 - 20.558591L) < 1e-6L);

    Real ans = S_value(50, depth_limit);
    std::cout << std::fixed << std::setprecision(6) << ans << '\n';
    return 0;
}
