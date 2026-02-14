#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

template <typename T>
struct GcdResult {
    T gcd;
    T coeff_a;
    T coeff_b;
};

template <typename T>
static GcdResult<T> extended_gcd(T a, T b) {
    T x = a, y = b;
    T ax = 1, ay = 0;
    T bx = 0, by = 1;
    while (x != 0) {
        const T k = y / x;
        y %= x;
        ay -= k * ax;
        by -= k * bx;
        std::swap(x, y);
        std::swap(ax, ay);
        std::swap(bx, by);
    }
    return {y, ay, by};
}

struct GaussianInt {
    int x;
    int y;

    friend GaussianInt operator*(const GaussianInt& a, const GaussianInt& b) {
        return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
    }

    int norm() const { return x * x + y * y; }
};

static int64_t solve_F(const int N) {
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(N / 2));
    std::vector<int> spf(static_cast<std::size_t>(N + 1), 0);
    for (int i = 2; i <= N; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (const int p : primes) {
            if (static_cast<int64_t>(p) * i > N) break;
            spf[static_cast<std::size_t>(p * i)] = p;
            if (p == spf[static_cast<std::size_t>(i)]) break;
        }
    }

    std::vector<int> exact_sqrt(static_cast<std::size_t>(N + 1), 0);
    for (int i = 1; static_cast<int64_t>(i) * i <= N; ++i) {
        exact_sqrt[static_cast<std::size_t>(i * i)] = i;
    }

    std::vector<std::optional<GaussianInt>> vals(static_cast<std::size_t>(N + 1));

    const auto normalize = [](GaussianInt a) -> GaussianInt {
        a.x = std::abs(a.x);
        a.y = std::abs(a.y);
        if (a.x < a.y) std::swap(a.x, a.y);
        assert(a.x >= a.y);
        return a;
    };

    vals[1] = GaussianInt{1, 0};
    if (N >= 2) vals[2] = GaussianInt{1, 1};

    for (const int p : primes) {
        if (p == 2 || (p & 3) == 3) continue;

        int a = 0, b = 0;
        for (int i = 1;; ++i) {
            const int rem = p - i * i;
            if (rem >= 0 && exact_sqrt[static_cast<std::size_t>(rem)] != 0) {
                a = i;
                b = exact_sqrt[static_cast<std::size_t>(rem)];
                break;
            }
        }
        assert(a <= b);

        GaussianInt cur{1, 0};
        for (int q = 1; static_cast<int64_t>(q) * p <= N;) {
            q *= p;
            cur = cur * GaussianInt{a, b};
            vals[static_cast<std::size_t>(q)] = normalize(cur);
            if (2 * q <= N) {
                vals[static_cast<std::size_t>(2 * q)] = normalize(cur * GaussianInt{1, 1});
            }
        }
    }

    std::vector<std::vector<int>> parents(static_cast<std::size_t>(N + 1));
    std::vector<std::vector<std::pair<int64_t, int64_t>>> dp(static_cast<std::size_t>(N + 1));
    for (int i = 1; i <= N; ++i) {
        if (!vals[static_cast<std::size_t>(i)].has_value()) continue;
        assert(vals[static_cast<std::size_t>(i)]->norm() == i);
        if (i == 1) continue;

        auto [g, cx, cy] = extended_gcd(vals[static_cast<std::size_t>(i)]->x, vals[static_cast<std::size_t>(i)]->y);
        if (g < 0) {
            g = -g;
            cx = -cx;
            cy = -cy;
        }
        assert(g == 1);

        GaussianInt p1{std::abs(cy), std::abs(cx)};
        GaussianInt p2{vals[static_cast<std::size_t>(i)]->x - p1.x, vals[static_cast<std::size_t>(i)]->y - p1.y};
        if (p1.norm() > p2.norm()) std::swap(p1, p2);
        if (p1.norm() == p2.norm()) assert(i == 2);

        if (vals[static_cast<std::size_t>(p1.norm())].has_value()) {
            parents[static_cast<std::size_t>(i)].push_back(p1.norm());
        }
        if (i > 2 && vals[static_cast<std::size_t>(p2.norm())].has_value()) {
            parents[static_cast<std::size_t>(i)].push_back(p2.norm());
        }

        dp[static_cast<std::size_t>(i)].assign(parents[static_cast<std::size_t>(i)].size(), {0, 0});
    }

    int64_t ways = 0;
    int64_t answer = 0;
    for (int i = N; i >= 1; --i) {
        auto& pi = parents[static_cast<std::size_t>(i)];
        auto& dpi = dp[static_cast<std::size_t>(i)];

        for (int z = 0; z < static_cast<int>(pi.size()); ++z) {
            auto& [cur_ways, cur_sum] = dpi[static_cast<std::size_t>(z)];
            ways += cur_ways;
            answer += cur_ways * static_cast<int64_t>(i + pi[static_cast<std::size_t>(z)]) + cur_sum;
            ++cur_ways;
        }

        assert(pi.size() <= 2);
        if (pi.size() == 2) {
            const int j = pi[1];
            auto& pj = parents[static_cast<std::size_t>(j)];
            const int z = static_cast<int>(
                std::find(pj.begin(), pj.end(), pi[0]) - pj.begin()
            );
            assert(z < static_cast<int>(pj.size()) && pj[static_cast<std::size_t>(z)] == pi[0]);

            const auto [w0, s0] = dpi[0];
            const auto [w1, s1] = dpi[1];
            auto& [wj, sj] = dp[static_cast<std::size_t>(j)][static_cast<std::size_t>(z)];
            wj += w0 * w1;
            sj += w0 * w1 * static_cast<int64_t>(i) + w0 * s1 + w1 * s0;
        }
    }

    (void)ways;
    return answer;
}

int main() {
    assert(solve_F(20) == 258);
    assert(solve_F(100) == 538768);
    std::cout << solve_F(1'000'000) << '\n';
    return 0;
}
