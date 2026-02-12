#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;

long double central_ratio_asymptotic(const u64 n) {
    // C(2n,n) / 4^n = 1/sqrt(pi*n) * (1 - 1/(8n) + 1/(128n^2) + 5/(1024n^3) - ...).
    const long double nn = static_cast<long double>(n);
    const long double inv = 1.0L / nn;
    const long double corr = 1.0L - inv / 8.0L
                           + (inv * inv) / 128.0L
                           + (5.0L * inv * inv * inv) / 1024.0L
                           - (21.0L * inv * inv * inv * inv) / 32768.0L;
    return corr / std::sqrtl(std::acos(-1.0L) * nn);
}

u64 g_from_threshold(const long double threshold) {
    // For small n, exact recurrence from r_0 = 1 and
    // r_n = r_{n-1} * (2n-1)/(2n).
    const long double pi = std::acos(-1.0L);
    const u64 approx = static_cast<u64>(std::floor(1.0L / (pi * threshold * threshold)));
    if (approx < 100'000ULL) {
        u64 n = 0;
        long double r = 1.0L;
        while (r > threshold) {
            ++n;
            r *= (2.0L * static_cast<long double>(n) - 1.0L) / (2.0L * static_cast<long double>(n));
        }
        return n;
    }

    u64 n = std::max<u64>(1ULL, approx);
    long double r = central_ratio_asymptotic(n);

    while (n > 1 && r <= threshold) {
        r *= (2.0L * static_cast<long double>(n)) / (2.0L * static_cast<long double>(n) - 1.0L);
        --n;
    }

    while (r > threshold) {
        ++n;
        r *= (2.0L * static_cast<long double>(n) - 1.0L) / (2.0L * static_cast<long double>(n));
    }

    return n;
}

u64 g_from_x(const long double x) {
    const long double threshold = 2.0L / x - 1.0L;
    return g_from_threshold(threshold);
}

}  // namespace

int main() {
    assert(g_from_x(1.7L) == 10ULL);

    // For X = 1.9999 = 19999/10000, threshold is exactly 1/19999.
    constexpr long double threshold = 1.0L / 19'999.0L;
    std::cout << g_from_threshold(threshold) << '\n';
    return 0;
}
