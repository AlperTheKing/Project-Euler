#include <boost/multiprecision/cpp_dec_float.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Dec = boost::multiprecision::cpp_dec_float_100;

std::uint64_t choose_int(int n, int r) {
    if (r < 0 || r > n) {
        return 0;
    }
    r = std::min(r, n - r);
    std::uint64_t result = 1;
    for (int i = 1; i <= r; ++i) {
        result = (result * static_cast<std::uint64_t>(n - r + i)) / static_cast<std::uint64_t>(i);
    }
    return result;
}

Dec factorial_dec(int n) {
    Dec result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= Dec(i);
    }
    return result;
}

Dec probability_cycle_shift(int k) {
    const int n = k * (k - 1) / 2 + 1;

    std::vector<std::vector<Dec>> a(k + 1, std::vector<Dec>(n, Dec(1)));

    for (int i = 1; i <= k; ++i) {
        std::vector<Dec> expected_fixed_subsets(n + 1, Dec(0));

        for (int r = 0; r <= n; ++r) {
            Dec p_invariant = 0;
            const int j_min = std::max(0, i - (n - r));
            const int j_max = std::min(i, r);

            for (int j = j_min; j <= j_max; ++j) {
                const Dec ways = Dec(choose_int(r, j)) * Dec(choose_int(n - r, i - j));
                const Dec prob_intersection_size = ways / Dec(choose_int(n, i));
                const Dec prob_preserve = Dec(1) / Dec(choose_int(i, j));
                p_invariant += prob_intersection_size * prob_preserve;
            }

            expected_fixed_subsets[r] = Dec(choose_int(n, r)) * p_invariant;
        }

        for (int r = 0; r < n; ++r) {
            const Dec expected_hook_character =
                expected_fixed_subsets[r] - (r == 0 ? Dec(0) : expected_fixed_subsets[r - 1]);
            const Dec hook_dimension = Dec(choose_int(n - 1, r));
            a[i][r] = expected_hook_character / hook_dimension;
        }
    }

    Dec total = 0;
    for (int r = 0; r < n; ++r) {
        Dec product = 1;
        for (int i = 1; i <= k; ++i) {
            product *= a[i][r];
        }
        const Dec term = Dec(choose_int(n - 1, r)) * product;
        total += (r % 2 == 0 ? term : -term);
    }

    return total / factorial_dec(n);
}

bool approx_equal(const Dec& a, const Dec& b, const Dec& eps) {
    Dec diff = a - b;
    if (diff < 0) {
        diff = -diff;
    }
    return diff <= eps;
}

bool run_checkpoints() {
    {
        const Dec p2 = probability_cycle_shift(2);
        if (!approx_equal(p2, Dec("0.5"), Dec("1e-30"))) {
            std::cerr << "Checkpoint failed for P(2): got " << p2 << '\n';
            return false;
        }
    }

    {
        const Dec p3 = probability_cycle_shift(3);
        const Dec expected = Dec(1) / Dec(72);
        if (!approx_equal(p3, expected, Dec("1e-30"))) {
            std::cerr << "Checkpoint failed for P(3): got " << p3 << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    const Dec answer = probability_cycle_shift(7);
    std::cout << std::scientific << std::setprecision(10) << answer << '\n';
    return 0;
}
