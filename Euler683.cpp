#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using ld = long double;
using cd = std::complex<ld>;

ld expected_round_payment(int m) {
    std::vector<ld> den(static_cast<std::size_t>(m), 0.0L);
    std::vector<ld> mu(static_cast<std::size_t>(m), 0.0L);

    for (int k = 1; k < m; ++k) {
        const ld th = 2.0L * acosl(-1.0L) * static_cast<ld>(k) / static_cast<ld>(m);
        const ld lam = 1.0L / 3.0L + (4.0L / 9.0L) * cosl(th) + (2.0L / 9.0L) * cosl(2.0L * th);
        den[static_cast<std::size_t>(k)] = 1.0L - lam;
    }

    for (int k = 1; k < m; ++k) {
        const ld inv = 1.0L / den[static_cast<std::size_t>(k)];
        const ld th = 2.0L * acosl(-1.0L) * static_cast<ld>(k) / static_cast<ld>(m);
        const cd base(cosl(th), sinl(th));
        cd cur(1.0L, 0.0L);

        for (int d = 1; d < m; ++d) {
            cur *= base;
            mu[static_cast<std::size_t>(d)] += (1.0L - cur.real()) * inv;
        }
    }

    std::vector<ld> r(static_cast<std::size_t>(m), 0.0L);
    ld sum_r = 0.0L;
    for (int d = 1; d < m; ++d) {
        r[static_cast<std::size_t>(d)] = 2.0L * mu[static_cast<std::size_t>(d)] - 1.0L;
        sum_r += r[static_cast<std::size_t>(d)];
    }
    r[0] = -sum_r;

    cd sum_s(0.0L, 0.0L);
    for (int k = 1; k < m; ++k) {
        const ld th = -2.0L * acosl(-1.0L) * static_cast<ld>(k) / static_cast<ld>(m);
        const cd base(cosl(th), sinl(th));
        cd cur(1.0L, 0.0L);
        cd rhat(0.0L, 0.0L);

        for (int d = 0; d < m; ++d) {
            rhat += cur * r[static_cast<std::size_t>(d)];
            cur *= base;
        }

        sum_s += rhat / den[static_cast<std::size_t>(k)];
    }

    return -sum_s.real() / static_cast<ld>(m);
}

std::vector<ld> G_prefix(int n) {
    std::vector<ld> pref(static_cast<std::size_t>(n + 1), 0.0L);
    ld acc = 0.0L;
    for (int m = 2; m <= n; ++m) {
        acc += expected_round_payment(m);
        pref[static_cast<std::size_t>(m)] = acc;
    }
    return pref;
}

ld round_sig(ld x, int sig) {
    const ld e = floorl(log10l(x));
    const ld scale = powl(10.0L, static_cast<ld>(sig - 1) - e);
    return roundl(x * scale) / scale;
}

}  // namespace

int main() {
    const std::vector<ld> g = G_prefix(500);

    assert(fabsl(g[5] - 96.544L) < 0.001L);
    assert(fabsl(round_sig(g[50], 9) - 2.82491788e6L) < 0.5L);

    std::cout << std::scientific << std::setprecision(8) << static_cast<double>(g[500]) << "\n";
    return 0;
}
