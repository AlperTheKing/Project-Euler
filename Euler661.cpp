#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {

long double expected_leads(long double pA, long double pB, long double pStop) {
    const long double q = 1.0L - pStop;
    const long double c = 1.0L - q * (1.0L - pA - pB);
    const long double disc = c * c - 4.0L * q * q * pA * pB;
    const long double root = std::sqrt(disc);

    const long double r_small = (c - root) / (2.0L * q * pA);
    const long double r_big = (c + root) / (2.0L * q * pA);

    return (1.0L - r_small) / (pStop * q * (r_big - r_small));
}

long double solve() {
    const long double e1 = expected_leads(0.25L, 0.25L, 0.5L);
    const long double e2 = expected_leads(0.47L, 0.48L, 0.001L);
    assert(std::fabsl(e1 - 0.585786L) < 1e-6L);
    assert(std::fabsl(e2 - 377.471736L) < 1e-6L);

    long double h3 = 0.0L;
    {
        const int k = 3;
        const long double pA = 1.0L / std::sqrt(static_cast<long double>(k + 3));
        const long double pB = pA + 1.0L / static_cast<long double>(k * k);
        const long double p = 1.0L / static_cast<long double>(k) / static_cast<long double>(k) /
                              static_cast<long double>(k);
        h3 = expected_leads(pA, pB, p);
    }
    assert(std::fabsl(h3 - 6.8345L) < 1e-4L);

    long double ans = 0.0L;
    for (int k = 3; k <= 50; ++k) {
        const long double pA = 1.0L / std::sqrt(static_cast<long double>(k + 3));
        const long double pB = pA + 1.0L / static_cast<long double>(k * k);
        const long double p = 1.0L / static_cast<long double>(k) / static_cast<long double>(k) /
                              static_cast<long double>(k);
        ans += expected_leads(pA, pB, p);
    }
    return ans;
}

}  // namespace

int main() {
    std::cout << std::fixed << std::setprecision(4) << solve() << "\n";
    return 0;
}
