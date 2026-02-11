#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

template <class F>
long double expected_error(const int n, const int m, F f_value) {
    if (m >= n) {
        return 0.0L;
    }

    long double p = static_cast<long double>(n - m) / static_cast<long double>(n);
    long double total = static_cast<long double>(f_value(1)) * p;

    for (int k = 1; k < n - m; ++k) {
        p *= static_cast<long double>(n - k - m) / static_cast<long double>(n - k);
        total += static_cast<long double>(f_value(k + 1)) * p;
    }

    return total;
}

std::vector<int> totients_up_to(const int n) {
    std::vector<int> phi(static_cast<std::size_t>(n + 1));
    for (int i = 0; i <= n; ++i) {
        phi[static_cast<std::size_t>(i)] = i;
    }
    for (int p = 2; p <= n; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (int j = p; j <= n; j += p) {
            phi[static_cast<std::size_t>(j)] -= phi[static_cast<std::size_t>(j)] / p;
        }
    }
    return phi;
}

}  // namespace

int main() {
    {
        const long double sample = expected_error(100, 50, [](const int k) { return k; });
        const long double exact = 2525.0L / 1326.0L;
        assert(std::fabsl(sample - exact) < 1e-12L);
    }

    {
        const std::vector<int> phi = totients_up_to(10'000);
        const long double sample =
            expected_error(10'000, 100, [&](const int k) { return phi[static_cast<std::size_t>(k)]; });
        assert(std::fabsl(sample - 5842.849907L) < 1e-6L);
    }

    const int n = 12'345'678;
    const int m = 12'345;
    const std::vector<int> phi = totients_up_to(n);
    const long double answer =
        expected_error(n, m, [&](const int k) { return phi[static_cast<std::size_t>(k)]; });

    std::cout << std::fixed << std::setprecision(6) << answer << '\n';
    return 0;
}
