#include <boost/multiprecision/cpp_dec_float.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>

using Dec = boost::multiprecision::cpp_dec_float_100;

static std::uint64_t choose_u64(int n, int r) {
    if (r < 0 || r > n) return 0;
    r = std::min(r, n - r);
    std::uint64_t x = 1;
    for (int i = 1; i <= r; ++i) {
        x = (x * static_cast<std::uint64_t>(n - r + i)) / static_cast<std::uint64_t>(i);
    }
    return x;
}

static Dec factorial_dec(int n) {
    Dec x = 1;
    for (int i = 2; i <= n; ++i) x *= Dec(i);
    return x;
}

static Dec probability(int k) {
    const int n = k * (k - 1) / 2 + 1;
    Dec sum = 0;
    for (int r = 0; r < n; ++r) {
        const Dec dim = Dec(choose_u64(n - 1, r));
        Dec prod = 1;
        for (int i = 1; i <= k; ++i) {
            prod *= Dec(choose_u64(n - i, r)) / dim;
        }
        Dec term = dim * prod;
        if (r & 1) term = -term;
        sum += term;
    }
    return sum / factorial_dec(n);
}

static bool close_to(const Dec& a, const Dec& b, const Dec& eps) {
    Dec d = a - b;
    if (d < 0) d = -d;
    return d <= eps;
}

int main() {
    assert(close_to(probability(2), Dec("5e-1"), Dec("1e-40")));
    assert(close_to(probability(3), Dec("1.38888888888888888888888888888888888889e-2"), Dec("1e-40")));
    assert(close_to(probability(4), Dec("6.61375661375661375661375661375661375661e-6"), Dec("1e-40")));

    const Dec ans = probability(7);
    std::cout << std::scientific << std::setprecision(10) << ans << '\n';
    return 0;
}
