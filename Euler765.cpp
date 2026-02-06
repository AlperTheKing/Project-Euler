#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

long double optimal_success_probability(int rounds, long double p_win, long double target) {
    if (rounds < 0) {
        throw std::invalid_argument("rounds must be non-negative");
    }
    if (!(p_win > 0.0L && p_win < 1.0L)) {
        throw std::invalid_argument("p_win must be in (0, 1)");
    }
    if (target <= 1.0L) {
        return 1.0L;
    }

    const long double budget_q = 1.0L / target;

    // q[w] = Q(W = w), where Q is the fair-coin measure (0.5 / 0.5).
    std::vector<long double> q(rounds + 1, 0.0L);
    q[0] = std::ldexp(1.0L, -rounds);
    for (int w = 0; w < rounds; ++w) {
        q[w + 1] = q[w] * static_cast<long double>(rounds - w) / static_cast<long double>(w + 1);
    }

    // p[w] = P(W = w), where P uses the true win probability p_win.
    // Using p[w] = q[w] * (2(1-p))^rounds * (p/(1-p))^w keeps numbers stable.
    std::vector<long double> p(rounds + 1, 0.0L);
    const long double base = std::powl(2.0L * (1.0L - p_win), rounds);
    const long double ratio = p_win / (1.0L - p_win);
    long double ratio_pow = 1.0L;
    for (int w = 0; w <= rounds; ++w) {
        p[w] = q[w] * base * ratio_pow;
        ratio_pow *= ratio;
    }

    // Neyman-Pearson on leaf likelihood ratios: include all largest w, then possibly
    // a fractional part of one boundary layer to use Q-budget exactly.
    long double used_q = 0.0L;
    long double success_p = 0.0L;
    int boundary_w = -1;

    for (int w = rounds; w >= 0; --w) {
        const long double candidate_q = used_q + q[w];
        if (candidate_q <= budget_q + 1e-30L) {
            used_q = candidate_q;
            success_p += p[w];
        } else {
            boundary_w = w;
            break;
        }
    }

    if (boundary_w == -1) {
        return 1.0L;
    }

    long double frac = (budget_q - used_q) / q[boundary_w];
    frac = std::clamp(frac, 0.0L, 1.0L);
    success_p += frac * p[boundary_w];

    return success_p;
}

void validate() {
    const long double tol = 1e-15L;

    // 1 round, target 2: all-in gives success exactly with one win.
    {
        long double got = optimal_success_probability(1, 0.6L, 2.0L);
        long double want = 0.6L;
        if (std::fabsl(got - want) > tol) {
            throw std::runtime_error("Validation failed: N=1, T=2");
        }
    }

    // 2 rounds, target 4: need two net doublings, probability p^2.
    {
        long double got = optimal_success_probability(2, 0.6L, 4.0L);
        long double want = 0.36L;
        if (std::fabsl(got - want) > tol) {
            throw std::runtime_error("Validation failed: N=2, T=4");
        }
    }

    // If target <= starting capital, success is guaranteed.
    {
        long double got = optimal_success_probability(1000, 0.6L, 1.0L);
        long double want = 1.0L;
        if (std::fabsl(got - want) > tol) {
            throw std::runtime_error("Validation failed: T<=1 guarantee");
        }
    }
}

}  // namespace

int main() {
    validate();

    const int rounds = 1000;
    const long double p_win = 0.6L;
    const long double target = 1.0e12L;

    long double ans = optimal_success_probability(rounds, p_win, target);

    std::cout << std::fixed << std::setprecision(10)
              << static_cast<double>(ans) << '\n';
    return 0;
}
