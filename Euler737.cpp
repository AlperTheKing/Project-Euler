#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

namespace {

using i64 = std::int64_t;

constexpr long double kPi = 3.141592653589793238462643383279502884L;

long double harmonic_approx(const long double n) {
    const long double n2 = n * n;
    return std::logl(n) + 0.5772156649L + 1.0L / (2.0L * n) - 1.0L / (12.0L * n2) +
           1.0L / (120.0L * n2 * n2) - 1.0L / (252.0L * n2 * n2 * n2);
}

long double tan_phi(const long double n, const long double hn) {
    return std::sqrtl(n / hn - 0.25L) / (n + 0.5L);
}

long double atan_poly(const long double x) {
    const long double x2 = x * x;
    return x - x * x2 / 3.0L + x * x2 * x2 / 5.0L;
}

long double integral_block(const long double k, const long double f_n, const long double f_n_der2) {
    return 2.0L * k * f_n + f_n_der2 * k * k * k / 3.0L;
}

i64 coins_needed_exact(const i64 loops) {
    const long double target = 2.0L * kPi * static_cast<long double>(loops);
    long double harmonic = 1.0L;
    long double radius = 1.0L;
    long double tan_gamma = 0.0L;
    long double accumulated = 0.0L;
    i64 n = 1;

    while (accumulated <= target) {
        const long double r2 = radius * radius;
        const long double root = std::sqrtl(4.0L - r2);
        const long double tan_beta = root / radius;
        const long double tan_theta = (tan_beta - tan_gamma) / (1.0L + tan_beta * tan_gamma);
        accumulated += std::atanl(tan_theta);

        const long double nl = static_cast<long double>(n);
        tan_gamma = (nl * radius * root) / (nl * r2 + 2.0L);

        const i64 m = n + 1;
        harmonic += 1.0L / static_cast<long double>(m);
        radius = std::sqrtl(harmonic / static_cast<long double>(m));
        ++n;
    }

    return n;
}

i64 coins_needed_approx(const i64 loops) {
    long double total_angle = -0.16103076705762L / 180.0L * kPi;
    i64 n = 0;
    long double harmonic = 0.0L;
    constexpr i64 initial = 1'000'000;
    constexpr i64 k = 20'000;

    for (i64 i = 0; i < initial; ++i) {
        ++n;
        harmonic += 1.0L / static_cast<long double>(n);
        total_angle += atan_poly(tan_phi(static_cast<long double>(n), harmonic));
    }

    i64 center = n + k + 1;
    long double block_sum = 0.0L;
    const long double stage1_target = static_cast<long double>(loops) * 2.0L * kPi - kPi / 2.0L;

    while (total_angle < stage1_target) {
        const long double f_nm2k =
            atan_poly(tan_phi(static_cast<long double>(center - 2 * k), harmonic_approx(center - 2.0L * k)));
        const long double f_nmk =
            atan_poly(tan_phi(static_cast<long double>(center - k), harmonic_approx(center - 1.0L * k)));
        const long double f_n = atan_poly(tan_phi(static_cast<long double>(center), harmonic_approx(center)));
        const long double f_npk =
            atan_poly(tan_phi(static_cast<long double>(center + k), harmonic_approx(center + 1.0L * k)));
        const long double f_np2k =
            atan_poly(tan_phi(static_cast<long double>(center + 2 * k), harmonic_approx(center + 2.0L * k)));

        const long double f_n_der2 = (f_npk + f_nmk - 2.0L * f_n) / (static_cast<long double>(k) * k);
        const long double f_npk_der = (-f_n + f_np2k) / (2.0L * static_cast<long double>(k));
        const long double f_nmk_der = (f_n - f_nm2k) / (2.0L * static_cast<long double>(k));
        block_sum = integral_block(static_cast<long double>(k), f_n, f_n_der2) + (f_nmk + f_npk) / 2.0L +
                    (f_npk_der - f_nmk_der) / 12.0L;

        total_angle += block_sum;
        center += 2 * k + 1;
    }

    center -= 2 * k + 1;
    total_angle -= block_sum;
    n = center - k;
    harmonic = harmonic_approx(static_cast<long double>(n));

    long double d_theta = std::atanl(std::sqrtl(4.0L * static_cast<long double>(n) / harmonic - 1.0L));
    long double d_phi = atan_poly(tan_phi(static_cast<long double>(n), harmonic));
    const long double stage2_target = static_cast<long double>(loops) * 2.0L * kPi;

    while (total_angle + d_theta < stage2_target) {
        total_angle += d_phi;
        ++n;
        harmonic += 1.0L / static_cast<long double>(n);
        d_phi = atan_poly(tan_phi(static_cast<long double>(n), harmonic));
        d_theta = std::atanl(std::sqrtl(4.0L * static_cast<long double>(n) / harmonic - 1.0L));
    }

    return n + 1;
}

i64 coins_needed(const i64 loops) {
    if (loops <= 50) {
        return coins_needed_exact(loops);
    }
    return coins_needed_approx(loops);
}

}  // namespace

int main() {
    assert(coins_needed(1) == 31);
    assert(coins_needed(2) == 154);
    assert(coins_needed(10) == 6947);
    std::cout << coins_needed(2020) << '\n';
    return 0;
}
