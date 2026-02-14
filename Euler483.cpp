#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

static inline int gcd_int(const int a, const int b) {
    return b == 0 ? a : gcd_int(b, a % b);
}

static inline int lcm_int(const int a, const int b) {
    return a / gcd_int(a, b) * b;
}

static double g_value(const int n) {
    std::vector<int> largest_prime_factor(static_cast<std::size_t>(n + 1), 0);
    for (int p = 2; p <= n; ++p) {
        if (largest_prime_factor[static_cast<std::size_t>(p)] != 0) continue;
        for (int v = p; v <= n; v += p) largest_prime_factor[static_cast<std::size_t>(v)] = p;
    }

    std::vector<std::map<int, double>> histograms(static_cast<std::size_t>(n + 1));
    {
        double count = 1.0;
        for (int i = 0; i <= n; ++i) {
            if (i > 0) count /= static_cast<double>(i);
            histograms[static_cast<std::size_t>(i)][1] = count;
        }
    }

    for (int p = n; p >= 2; --p) {
        if (largest_prime_factor[static_cast<std::size_t>(p)] != p) continue;

        for (int c = p; c <= n; c += p) {
            if (largest_prime_factor[static_cast<std::size_t>(c)] != p) continue;

            std::vector<std::map<int, double>> next(static_cast<std::size_t>(n + 1));
            double count_rhs = 1.0;
            for (int multiplicity = 0; multiplicity * c <= n;
                 ++multiplicity, count_rhs /= static_cast<double>(c) * static_cast<double>(multiplicity)) {
                const int rhs_len = multiplicity * c;
                const int rhs_lcm = (multiplicity == 0 ? 1 : c);

                for (int lhs_len = 0; lhs_len + rhs_len <= n; ++lhs_len) {
                    const int total = lhs_len + rhs_len;
                    const auto& h = histograms[static_cast<std::size_t>(lhs_len)];
                    auto& nh = next[static_cast<std::size_t>(total)];

                    for (const auto& kv : h) {
                        const int l = lcm_int(kv.first, rhs_lcm);
                        nh[l] += kv.second * count_rhs;
                    }
                }
            }
            histograms.swap(next);
        }

        for (int len = 0; len <= n; ++len) {
            auto& h = histograms[static_cast<std::size_t>(len)];
            std::map<int, double> squashed;
            for (const auto& kv : h) {
                int l = kv.first;
                double c = kv.second;
                while (l % p == 0) {
                    l /= p;
                    c *= static_cast<double>(p) * static_cast<double>(p);
                }
                squashed[l] += c;
            }
            h.swap(squashed);
        }
    }

    double ans = 0.0;
    for (const auto& kv : histograms[static_cast<std::size_t>(n)]) {
        ans += static_cast<double>(kv.first) * static_cast<double>(kv.first) * kv.second;
    }
    return ans;
}

static bool close_to(const double a, const double b, const double eps) {
    return std::fabs(a - b) <= eps;
}

int main() {
    assert(close_to(g_value(3), 31.0 / 6.0, 1e-12));
    assert(close_to(g_value(5), 2081.0 / 120.0, 1e-12));
    assert(close_to(g_value(20), 5106.136147, 2e-6));

    const double ans = g_value(350);
    std::cout << std::scientific << std::setprecision(9) << ans << '\n';
    return 0;
}
