#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

bool is_power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

bool is_one_less_than_power_of_two(int n) {
    return n >= 1 && ((n + 1) & n) == 0;
}

int highest_power_of_two_at_most(int n) {
    int p = 1;
    while ((p << 1) <= n) {
        p <<= 1;
    }
    return p;
}

long double log_u64(u64 x) {
    return std::log(static_cast<long double>(x));
}

long double log_k_infty(int n) {
    if (n == 0) {
        return (log_u64(2) + 2 * log_u64(4) + log_u64(6)) / 4;
    }

    const u64 b = (u64{1} << n) - 1;

    if (is_one_less_than_power_of_two(n)) {
        return (log_u64(2) + 2 * log_u64(b)) / 5;
    }

    if (is_power_of_two(n)) {
        return log_u64(b) / 2 + (log_u64(b - 1) + log_u64(b + 1)) / 12;
    }

    const int p = highest_power_of_two_at_most(n);
    const int q = p << 1;
    const int e = q - n;
    const u64 c = (u64{1} << e) - 1;

    return log_u64(b) / 3 + log_u64(c) / 6 + (log_u64(c - 1) + log_u64(c + 1)) / 12;
}

long double empirical_log_k_from_truncation(int n, int m) {
    const int M = 1 << m;

    cpp_int p = 0;
    for (int i = 0; i <= m; ++i) {
        p += cpp_int(1) << (M - (1 << i));
    }

    cpp_int q = cpp_int(1) << (M - n);
    p %= q;

    long double sum_logs = 0;
    std::uint64_t len = 0;

    while (p != 0) {
        const cpp_int a = q / p;
        const long double av = a.convert_to<long double>();
        sum_logs += std::log(av);
        ++len;

        const cpp_int r = q - a * p;
        q = p;
        p = r;
    }

    return sum_logs / static_cast<long double>(len);
}

void validate() {
    assert(std::fabs(std::exp(log_k_infty(2)) - 2.059767L) < 1e-6L);

    struct Check {
        int n;
        int m;
        long double rel_tol;
    };

    const std::array<Check, 8> checks = {{
        {0, 16, 6e-4L},
        {1, 16, 5e-5L},
        {2, 16, 5e-5L},
        {5, 16, 5e-5L},
        {20, 16, 3e-4L},
        {31, 16, 3e-3L},
        {32, 16, 3e-3L},
        {50, 16, 2e-3L},
    }};

    for (const auto& c : checks) {
        const long double expected = log_k_infty(c.n);
        const long double observed = empirical_log_k_from_truncation(c.n, c.m);
        const long double rel = std::fabs(std::exp(observed - expected) - 1.0L);
        assert(rel < c.rel_tol);
    }
}

}  // namespace

int main() {
    validate();

    long double sum_logs = 0;
    for (int n = 0; n <= 50; ++n) {
        sum_logs += log_k_infty(n);
    }

    const long double answer = std::exp(sum_logs / 51.0L);
    std::cout << std::fixed << std::setprecision(6) << answer << '\n';
    return 0;
}
