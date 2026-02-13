#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using ld = long double;

ld F(int n, int m) {
    std::vector<std::vector<ld>> binom(m + 1, std::vector<ld>(m + 1, 0.0L));
    for (int i = 0; i <= m; ++i) {
        binom[i][0] = 1.0L;
        binom[i][i] = 1.0L;
        for (int j = 1; j < i; ++j) {
            binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
        }
    }

    const ld two_pi = 2.0L * std::acos(-1.0L);
    std::vector<ld> cosine(n, 0.0L);
    for (int r = 0; r < n; ++r) {
        cosine[r] = std::cos(two_pi * static_cast<ld>(r) / static_cast<ld>(n));
    }

    ld total = 0.0L;
    auto dfs = [&](auto&& self, int idx, int rem, int mod_sum, int count_zero, ld cos_sum, ld ways) -> void {
        if (idx == n - 1) {
            const int x = rem;
            if ((mod_sum + idx * x) % n != 0) {
                return;
            }
            if (count_zero == m) {
                return;
            }
            const ld final_cos = cos_sum + static_cast<ld>(x) * cosine[idx];
            const ld denom = 1.0L - final_cos / static_cast<ld>(m);
            assert(denom > 0.0L);
            total += ways / denom;
            return;
        }

        for (int x = 0; x <= rem; ++x) {
            int next_count_zero = count_zero;
            if (idx == 0) {
                next_count_zero = x;
            }
            self(self, idx + 1, rem - x, (mod_sum + idx * x) % n, next_count_zero,
                 cos_sum + static_cast<ld>(x) * cosine[idx], ways * binom[rem][x]);
        }
    };

    dfs(dfs, 0, m, 0, 0, 0.0L, 1.0L);
    return total;
}

ld G(int N, int M) {
    ld total = 0.0L;
    for (int n = 2; n <= N; ++n) {
        for (int m = 2; m <= M; ++m) {
            total += F(n, m);
        }
    }
    return total;
}

void run_validations() {
    auto check = [](ld got, ld expected, ld tol) {
        assert(std::fabsl(got - expected) <= tol);
    };

    check(F(2, 2), 0.5L, 1e-15L);
    check(F(3, 2), 4.0L / 3.0L, 1e-15L);
    check(F(2, 3), 9.0L / 4.0L, 1e-15L);
    check(F(4, 5), 6875.0L / 24.0L, 1e-12L);
    check(G(3, 3), 137.0L / 12.0L, 1e-12L);
    check(G(4, 5), 6277.0L / 12.0L, 1e-10L);
    check(G(6, 6), 1.681521567954e4L, 1e-6L);
}

}  // namespace

int main() {
    run_validations();
    const ld answer = G(12, 12);
    std::cout << std::scientific << std::setprecision(12) << answer << '\n';
    return 0;
}
