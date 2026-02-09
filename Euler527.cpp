#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

long double harmonic_asymptotic(const std::uint64_t n) {
    // H_n = log(n) + gamma + 1/(2n) - 1/(12n^2) + 1/(120n^4) - ...
    // For n=1e10, terms past 1/n are negligible for 1e-8 output accuracy.
    static constexpr long double kGamma = 0.577215664901532860606512090082402431L;
    const long double nn = static_cast<long double>(n);
    const long double inv = 1.0L / nn;
    const long double inv2 = inv * inv;
    const long double inv4 = inv2 * inv2;
    return std::log(nn) + kGamma + 0.5L * inv - (1.0L / 12.0L) * inv2 + (1.0L / 120.0L) * inv4;
}

long double harmonic_number(const std::uint64_t n) {
    // Exact for small n; asymptotic for large n.
    if (n <= 1'000'000ULL) {
        long double s = 0.0L;
        for (std::uint64_t k = 1; k <= n; ++k) {
            s += 1.0L / static_cast<long double>(k);
        }
        return s;
    }
    return harmonic_asymptotic(n);
}

long double expected_random_binary_search(const std::uint64_t n) {
    // Random pivot each step is equivalent to searching for a random key in a random BST:
    // R(n) = 2*(n+1)/n * H_n - 3.
    const long double Hn = harmonic_number(n);
    const long double nn = static_cast<long double>(n);
    return 2.0L * (nn + 1.0L) / nn * Hn - 3.0L;
}

long double expected_standard_binary_search(const std::uint64_t n) {
    // Deterministic midpoint recursion gives a complete binary tree shape.
    // With h = floor(log2 n) and r = n - (2^h - 1), total comparisons are:
    //   sum_{d=0..h-1} (d+1)2^d + (h+1)r = (h+1)n - 2^{h+1} + (h+2).
    // B(n) is total / n.
    const int h = 63 - __builtin_clzll(n);  // floor(log2 n), n>0
    const std::uint64_t two_h1 = 1ULL << (h + 1);
    const long double nn = static_cast<long double>(n);
    return static_cast<long double>(h + 1) -
           (static_cast<long double>(two_h1) - static_cast<long double>(h + 2)) / nn;
}

bool run_checkpoints() {
    // Given examples.
    {
        const std::uint64_t n = 6;
        const long double B = expected_standard_binary_search(n);
        const long double R = expected_random_binary_search(n);
        if (std::fabsl(B - 2.3333333333333333L) > 1e-12L) {
            std::cerr << "Checkpoint failed: B(6)\n";
            return false;
        }
        if (std::fabsl(R - 2.7166666666666667L) > 1e-12L) {
            std::cerr << "Checkpoint failed: R(6)\n";
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

    constexpr std::uint64_t n = 10'000'000'000ULL;
    const long double B = expected_standard_binary_search(n);
    const long double R = expected_random_binary_search(n);
    const long double ans = R - B;
    std::cout << std::fixed << std::setprecision(8) << ans << '\n';
    return 0;
}
