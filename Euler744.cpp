#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

long double f_exact_small(const std::uint64_t n, const long double p) {
    const long double q = 1.0L - p;
    long double tv = std::powl(q, static_cast<long double>(n));
    long double te = std::powl(p, static_cast<long double>(n));

    long double total = 0.0L;
    for (std::uint64_t m = 0; m < n; ++m) {
        const long double weight = static_cast<long double>(n + 1ULL - m);
        total += weight * (tv + te);

        const long double mult = static_cast<long double>(n + m) /
                                 static_cast<long double>(m + 1ULL);
        tv *= mult * p;
        te *= mult * q;
    }

    return total / static_cast<long double>(2ULL * n + 1ULL);
}

long double f_large_drift(const std::uint64_t n, const long double p) {
    const long double q = 1.0L - p;
    const long double dominant = (p > q) ? p : q;

    const long double two_n_plus_one = static_cast<long double>(2ULL * n + 1ULL);
    const long double expected_questions = static_cast<long double>(n) / dominant;
    return (two_n_plus_one - expected_questions) / two_n_plus_one;
}

long double f(const std::uint64_t n, const long double p) {
    const long double q = 1.0L - p;
    const long double drift_score =
        std::fabsl(1.0L - 2.0L * p) *
        std::sqrt(static_cast<long double>(n) / (p * q));

    if (n <= 2000ULL || drift_score < 15.0L) {
        return f_exact_small(n, p);
    }
    return f_large_drift(n, p);
}

}  // namespace

int main() {
    assert(std::fabsl(f(6ULL, 0.5L) - 0.2851562500L) < 5e-11L);
    assert(std::fabsl(f(10ULL, 3.0L / 7.0L) - 0.2330040743L) < 5e-11L);
    assert(std::fabsl(f(10'000ULL, 0.3L) - 0.2857499982L) < 5e-10L);

    std::cout << std::fixed << std::setprecision(10)
              << static_cast<double>(f(100'000'000'000ULL, 0.4999L)) << '\n';
    return 0;
}
