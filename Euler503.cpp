#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

long double solve_expected_score(const int n) {
    // V_t = minimal expected final score from round t onward (before observing round-t rank).
    // Base: V_n = average rank in final forced round = (n + 1) / 2.
    long double next_value = (static_cast<long double>(n) + 1.0L) / 2.0L;

    for (int t = n - 1; t >= 1; --t) {
        const long double slope =
            (static_cast<long double>(n) + 1.0L) / (static_cast<long double>(t) + 1.0L);
        const long double cutoff = next_value / slope;

        int k = static_cast<int>(std::floor(cutoff + 1e-18L));
        if (k < 0) {
            k = 0;
        }
        if (k > t) {
            k = t;
        }

        const long double sum_ranks = static_cast<long double>(k) * (k + 1.0L) / 2.0L;
        const long double stop_expectation = slope * sum_ranks;
        const long double continue_expectation =
            static_cast<long double>(t - k) * next_value;

        next_value =
            (stop_expectation + continue_expectation) / static_cast<long double>(t);
    }
    return next_value;
}

bool nearly_equal(const long double a, const long double b, const long double eps) {
    return std::fabsl(a - b) <= eps;
}

bool run_checkpoints() {
    if (!nearly_equal(solve_expected_score(3), 5.0L / 3.0L, 1e-15L)) {
        std::cerr << "Checkpoint failed: F(3)\n";
        return false;
    }
    if (!nearly_equal(solve_expected_score(4), 15.0L / 8.0L, 1e-15L)) {
        std::cerr << "Checkpoint failed: F(4)\n";
        return false;
    }
    if (!nearly_equal(solve_expected_score(10), 2.5579365079L, 5e-11L)) {
        std::cerr << "Checkpoint failed: F(10)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr int n = 1'000'000;
    const long double answer = solve_expected_score(n);
    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(answer) << '\n';
    return 0;
}
