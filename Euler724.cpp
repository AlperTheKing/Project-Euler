#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

using i64 = std::int64_t;

constexpr long double kEulerGamma = 0.57721566490153286060651209008240243104215933593992L;

long double E_direct(const int n) {
    long double h1 = 0.0L;
    long double h2 = 0.0L;
    for (int i = 1; i <= n; ++i) {
        const long double inv = 1.0L / static_cast<long double>(i);
        h1 += inv;
        h2 += inv * inv;
    }
    return 0.5L * static_cast<long double>(n) * (h1 * h1 + h2);
}

long double E_large(const i64 n) {
    const long double nl = static_cast<long double>(n);
    const long double inv = 1.0L / nl;
    const long double inv2 = inv * inv;
    const long double inv3 = inv2 * inv;
    const long double inv4 = inv2 * inv2;
    const long double inv5 = inv4 * inv;
    const long double inv6 = inv3 * inv3;
    const long double inv7 = inv6 * inv;
    const long double inv8 = inv4 * inv4;
    const long double inv9 = inv8 * inv;

    const long double h1 = std::log(nl) + kEulerGamma +
                           0.5L * inv -
                           (1.0L / 12.0L) * inv2 +
                           (1.0L / 120.0L) * inv4 -
                           (1.0L / 252.0L) * inv6 +
                           (1.0L / 240.0L) * inv8;

    const long double pi = std::acos(-1.0L);
    const long double h2 = (pi * pi) / 6.0L -
                           inv +
                           0.5L * inv2 -
                           (1.0L / 6.0L) * inv3 +
                           (1.0L / 30.0L) * inv5 -
                           (1.0L / 42.0L) * inv7 +
                           (1.0L / 30.0L) * inv9;

    return 0.5L * nl * (h1 * h1 + h2);
}

}  // namespace

int main() {
    assert(std::fabsl(E_direct(2) - 3.5L) < 1e-12L);
    assert(std::fabsl(E_direct(5) - (12019.0L / 720.0L)) < 1e-12L);
    assert(std::fabsl(E_direct(100) - 1427.19347053489L) < 1e-9L);

    const long double ans = E_large(100'000'000LL);
    std::cout << std::llround(ans) << '\n';
    return 0;
}
