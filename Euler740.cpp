#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

struct Solver {
    int n;
    int dim;
    std::vector<long double> memo;
    std::vector<std::uint8_t> seen;

    explicit Solver(const int n_) : n(n_), dim(n_ + 1) {
        const std::size_t total = static_cast<std::size_t>(n + 1) * 3U * 3U *
                                  static_cast<std::size_t>(dim) * static_cast<std::size_t>(dim);
        memo.assign(total, 0.0L);
        seen.assign(total, 0U);
    }

    std::size_t index(const int i, const int x, const int y, const int a1, const int a2) const {
        std::size_t z = static_cast<std::size_t>(i);
        z = z * 3U + static_cast<std::size_t>(x);
        z = z * 3U + static_cast<std::size_t>(y);
        z = z * static_cast<std::size_t>(dim) + static_cast<std::size_t>(a1);
        z = z * static_cast<std::size_t>(dim) + static_cast<std::size_t>(a2);
        return z;
    }

    static void apply_draw(const int cat, int& g, int& y, int& a1, int& a2) {
        if (cat == 0) {
            --g;
        } else if (cat == 1) {
            --y;
        } else if (cat == 2) {
            --a1;
        } else {
            --a2;
            ++a1;
        }
    }

    long double solve(const int i, const int x, const int y, const int a1, const int a2) {
        if (i == n) {
            return (y > 0) ? 1.0L : 0.0L;
        }

        const std::size_t id = index(i, x, y, a1, a2);
        if (seen[id] != 0U) {
            return memo[id];
        }
        seen[id] = 1U;

        const int future_pool = n - i - 1;  // labels i+1..n-1
        const int total_slips = 2 * (n - i + 1);
        const int g0 = total_slips - x - y - a1 - 2 * a2;
        const int avail0 = total_slips - x;

        const auto count_for_cat = [](const int cat, const int g, const int yv, const int c1, const int c2) {
            if (cat == 0) return g;
            if (cat == 1) return yv;
            if (cat == 2) return c1;
            return 2 * c2;
        };

        long double ans = 0.0L;

        for (int c1 = 0; c1 < 4; ++c1) {
            const int cnt1 = count_for_cat(c1, g0, y, a1, a2);
            if (cnt1 <= 0) {
                continue;
            }

            const long double p1 = static_cast<long double>(cnt1) / static_cast<long double>(avail0);
            int g1 = g0, y1 = y, a11 = a1, a21 = a2;
            apply_draw(c1, g1, y1, a11, a21);

            const int avail1 = avail0 - 1;
            for (int c2 = 0; c2 < 4; ++c2) {
                const int cnt2 = count_for_cat(c2, g1, y1, a11, a21);
                if (cnt2 <= 0) {
                    continue;
                }

                const long double p2 = static_cast<long double>(cnt2) / static_cast<long double>(avail1);
                int g2 = g1, y2 = y1, a12 = a11, a22 = a21;
                apply_draw(c2, g2, y2, a12, a22);

                long double cont = 0.0L;
                if (i == n - 1) {
                    cont = (y2 > 0) ? 1.0L : 0.0L;
                } else {
                    const int b1 = a12;
                    const int b2 = a22;
                    const int b0 = future_pool - b1 - b2;

                    if (b0 > 0) {
                        const long double ps = static_cast<long double>(b0) / static_cast<long double>(future_pool);
                        cont += ps * solve(i + 1, 0, y2, b1, b2);
                    }
                    if (b1 > 0) {
                        const long double ps = static_cast<long double>(b1) / static_cast<long double>(future_pool);
                        cont += ps * solve(i + 1, 1, y2, b1 - 1, b2);
                    }
                    if (b2 > 0) {
                        const long double ps = static_cast<long double>(b2) / static_cast<long double>(future_pool);
                        cont += ps * solve(i + 1, 2, y2, b1, b2 - 1);
                    }
                }

                ans += p1 * p2 * cont;
            }
        }

        memo[id] = ans;
        return ans;
    }
};

long double q(const int n) {
    Solver solver(n);
    return solver.solve(1, 2, 2, 0, n - 2);
}

}  // namespace

int main() {
    assert(std::fabsl(q(3) - 0.3611111111L) < 5e-11L);
    assert(std::fabsl(q(5) - 0.2476095994L) < 5e-11L);

    std::cout << std::fixed << std::setprecision(10)
              << static_cast<double>(q(100)) << '\n';
    return 0;
}
