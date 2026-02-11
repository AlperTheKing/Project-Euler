#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

namespace {

using i64 = std::int64_t;

constexpr long double kPi = 3.141592653589793238462643383279502884L;

i64 coins_needed(const i64 loops) {
    const long double target = 2.0L * kPi * static_cast<long double>(loops);

    // After m coins are fixed, let G_m be their centroid in the table projection,
    // r_m = |G_m|, and gamma_m be the angle between the newest coin center and G_m.
    // The balancing geometry gives:
    // tan(theta_{m+1}) = (tan(beta_m)-tan(gamma_m)) / (1+tan(beta_m)tan(gamma_m)),
    // with cos(beta_m)=r_m/2 and
    // tan(gamma_{m+1}) = m*r_m*sqrt(4-r_m^2)/(m*r_m^2+2).
    // Also r_m^2 = H_m / m where H_m is harmonic sum.

    long double harmonic = 1.0L;   // H_1
    long double radius = 1.0L;     // r_1
    long double tan_gamma = 0.0L;  // gamma_1

    long double accumulated = 0.0L;
    i64 n = 1;

    while (accumulated <= target) {
        const long double root = std::sqrt(4.0L - radius * radius);
        const long double tan_beta = root / radius;
        const long double tan_theta =
            (tan_beta - tan_gamma) / (1.0L + tan_beta * tan_gamma);
        accumulated += std::atan(tan_theta);

        tan_gamma = (static_cast<long double>(n) * radius * root) /
                    (static_cast<long double>(n) * radius * radius + 2.0L);

        const i64 m = n + 1;
        harmonic += 1.0L / static_cast<long double>(m);
        radius = std::sqrt(harmonic / static_cast<long double>(m));
        ++n;
    }

    return n;
}

}  // namespace

int main() {
    assert(coins_needed(1) == 31);
    assert(coins_needed(2) == 154);
    assert(coins_needed(10) == 6947);

    std::cout << coins_needed(2020) << '\n';
    return 0;
}
