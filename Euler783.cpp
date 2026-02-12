#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

using ld = long double;

static ld solve_expectation(int n, int k) {
    ld m1 = static_cast<ld>(k);
    ld m2 = static_cast<ld>(k) * static_cast<ld>(k);
    ld total = 0;

    for (int t = 1; t <= n; ++t) {
        const ld m = static_cast<ld>(n - t + 2);
        const ld a = 2.0L / m;
        const ld c = (2.0L * static_cast<ld>(2 * k - 1)) /
                     (m * (static_cast<ld>(k) * m - 1.0L));

        total += c * m2 + (a - c) * m1;

        if (t == n) {
            break;
        }

        const ld r = (m - 2.0L) / m;
        const ld s = ((m - 2.0L) * (static_cast<ld>(k) * (m - 2.0L) - 1.0L)) /
                     (m * (static_cast<ld>(k) * m - 1.0L));

        const ld next_m2 = s * m2 + (2.0L * static_cast<ld>(k) * r + r - s) * m1 +
                           static_cast<ld>(k) * static_cast<ld>(k);
        const ld next_m1 = r * m1 + static_cast<ld>(k);

        m1 = next_m1;
        m2 = next_m2;
    }

    return total;
}

static ld log_binom(int n, int r) {
    if (r < 0 || r > n) {
        return -INFINITY;
    }
    return lgammal(static_cast<ld>(n) + 1.0L) -
           lgammal(static_cast<ld>(r) + 1.0L) -
           lgammal(static_cast<ld>(n - r) + 1.0L);
}

static ld hypergeom_prob(int population, int success, int draw, int picked_success) {
    if (picked_success < 0 || picked_success > success ||
        picked_success > draw || draw - picked_success > population - success) {
        return 0.0L;
    }
    const ld lp = log_binom(success, picked_success) +
                  log_binom(population - success, draw - picked_success) -
                  log_binom(population, draw);
    return expl(lp);
}

static ld brute_expectation_small(int n, int k) {
    std::vector<ld> dist(k * n + 1, 0.0L), next_dist(k * n + 1, 0.0L);
    dist[0] = 1.0L;

    ld answer = 0.0L;

    for (int t = 1; t <= n; ++t) {
        std::fill(next_dist.begin(), next_dist.end(), 0.0L);
        const int population = k * (n - t + 2);

        for (int y_prev = 0; y_prev < static_cast<int>(dist.size()); ++y_prev) {
            const ld prob_state = dist[y_prev];
            if (prob_state == 0.0L) {
                continue;
            }

            const int x = y_prev + k;
            for (int b = 0; b <= 2 * k; ++b) {
                const ld p = hypergeom_prob(population, x, 2 * k, b);
                if (p == 0.0L) {
                    continue;
                }
                answer += prob_state * p * static_cast<ld>(b) * static_cast<ld>(b);
                const int y_next = x - b;
                next_dist[y_next] += prob_state * p;
            }
        }

        dist.swap(next_dist);
    }

    return answer;
}

int main() {
    const ld sample = solve_expectation(2, 2);
    assert(fabsl(sample - 9.6L) < 1e-12L);

    const ld check_1 = solve_expectation(3, 2);
    const ld brute_1 = brute_expectation_small(3, 2);
    assert(fabsl(check_1 - brute_1) < 1e-11L);

    const ld check_2 = solve_expectation(5, 1);
    const ld brute_2 = brute_expectation_small(5, 1);
    assert(fabsl(check_2 - brute_2) < 1e-11L);

    const ld ans = solve_expectation(1'000'000, 10);
    std::cout << static_cast<long long>(llround(ans)) << '\n';
    return 0;
}
