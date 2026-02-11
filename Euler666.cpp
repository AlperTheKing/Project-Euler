#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

long double extinction_probability(int k, int m) {
    const int n = k * m;
    std::vector<int> r(static_cast<std::size_t>(n), 0);
    r[0] = 306;
    for (int i = 1; i < n; ++i) {
        const std::int64_t x = static_cast<std::int64_t>(r[static_cast<std::size_t>(i - 1)]);
        r[static_cast<std::size_t>(i)] = static_cast<int>((x * x) % 10'007LL);
    }

    std::vector<std::array<int, 5>> cnt(static_cast<std::size_t>(k));
    for (int i = 0; i < k; ++i) {
        cnt[static_cast<std::size_t>(i)].fill(0);
        for (int j = 0; j < m; ++j) {
            const int q = r[static_cast<std::size_t>(i * m + j)] % 5;
            ++cnt[static_cast<std::size_t>(i)][static_cast<std::size_t>(q)];
        }
    }

    std::vector<long double> p(static_cast<std::size_t>(k), 0.0L);
    std::vector<long double> next(static_cast<std::size_t>(k), 0.0L);

    constexpr int kMaxIter = 1'000'000;
    constexpr long double kEps = 1e-20L;

    for (int iter = 0; iter < kMaxIter; ++iter) {
        long double max_delta = 0.0L;

        for (int i = 0; i < k; ++i) {
            const auto& c = cnt[static_cast<std::size_t>(i)];
            const int a = (2 * i) % k;
            const int b = (i * i + 1) % k;
            const int d = (i + 1) % k;

            const long double pi = p[static_cast<std::size_t>(i)];
            const long double pa = p[static_cast<std::size_t>(a)];
            const long double pb = p[static_cast<std::size_t>(b)];
            const long double pd = p[static_cast<std::size_t>(d)];

            const long double value =
                (static_cast<long double>(c[0]) +
                 static_cast<long double>(c[1]) * pi * pi +
                 static_cast<long double>(c[2]) * pa +
                 static_cast<long double>(c[3]) * pb * pb * pb +
                 static_cast<long double>(c[4]) * pi * pd) /
                static_cast<long double>(m);

            next[static_cast<std::size_t>(i)] = value;
            const long double delta = std::fabsl(value - pi);
            if (delta > max_delta) {
                max_delta = delta;
            }
        }

        p.swap(next);
        if (max_delta < kEps) {
            break;
        }
    }

    return p[0];
}

std::int64_t rounded_1e8(long double x) {
    return static_cast<std::int64_t>(std::llround(x * 100'000'000.0L));
}

}  // namespace

int main() {
    assert(rounded_1e8(extinction_probability(2, 2)) == 7'243'802LL);
    assert(rounded_1e8(extinction_probability(4, 3)) == 18'554'021LL);
    assert(rounded_1e8(extinction_probability(10, 5)) == 53'466'253LL);

    const long double ans = extinction_probability(500, 10);
    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(ans) << "\n";
    return 0;
}
