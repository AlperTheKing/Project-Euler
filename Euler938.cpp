#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

long double solve(int R, int B) {
    std::vector<long double> prev(R + 1, 0.0L);
    std::vector<long double> curr(R + 1, 0.0L);

    for (int b = 1; b <= B; ++b) {
        curr[0] = 1.0L;
        for (int r = 1; r <= R; ++r) {
            if (r == 1) {
                curr[r] = prev[r];
                continue;
            }
            const long double numerator =
                static_cast<long double>(r - 1) * curr[r - 2] +
                static_cast<long double>(2 * b) * prev[r];
            const long double denominator = static_cast<long double>(r - 1 + 2 * b);
            curr[r] = numerator / denominator;
        }
        prev.swap(curr);
    }

    return prev[R];
}

void run_validations() {
    assert(std::fabsl(solve(2, 2) - 0.4666666667L) < 1e-9L);
    assert(std::fabsl(solve(10, 9) - 0.4118903397L) < 1e-9L);
    assert(std::fabsl(solve(34, 25) - 0.3665688069L) < 1e-9L);
}

}  // namespace

int main() {
    run_validations();
    std::cout << std::fixed << std::setprecision(10) << solve(24'690, 12'345) << '\n';
    return 0;
}
