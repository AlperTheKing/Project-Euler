#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;

static constexpr int MOD = 1'000'000'007;

static std::vector<int> distinct_parts(int limit) {
    std::vector<int> q(static_cast<std::size_t>(limit + 1), 0);
    q[0] = 1;

    std::vector<std::int8_t> pent_sign(static_cast<std::size_t>(limit + 1), 0);
    const int lim = static_cast<int>(std::sqrt(static_cast<double>(limit)));
    for (int j = 0; j <= lim; ++j) {
        const int sign = (j & 1) ? -1 : 1;
        const i64 g1 = static_cast<i64>(j) * (3LL * j + 1) / 2;
        const i64 g2 = static_cast<i64>(j) * (3LL * j - 1) / 2;
        if (g1 <= limit) pent_sign[static_cast<std::size_t>(g1)] = static_cast<std::int8_t>(sign);
        if (g2 <= limit) pent_sign[static_cast<std::size_t>(g2)] = static_cast<std::int8_t>(sign);
    }

    std::vector<int> squares;
    squares.reserve(static_cast<std::size_t>(lim));
    for (int k = 1; k <= lim; ++k) squares.push_back(k * k);

    for (int n = 1; n <= limit; ++n) {
        int val = pent_sign[static_cast<std::size_t>(n)];
        for (int k = 0; k < static_cast<int>(squares.size()) && squares[static_cast<std::size_t>(k)] <= n; ++k) {
            int term = q[static_cast<std::size_t>(n - squares[static_cast<std::size_t>(k)])];
            term <<= 1;
            if (term >= MOD) term -= MOD;
            if ((k & 1) == 0) {
                val += term;
                if (val >= MOD) val -= MOD;
            } else {
                val -= term;
                if (val < 0) val += MOD;
            }
        }
        q[static_cast<std::size_t>(n)] = val;
    }
    return q;
}

static std::vector<int> self_conjugate(int limit) {
    const int kmax = static_cast<int>(std::sqrt(static_cast<double>(limit)));
    std::vector<int> q_odd(static_cast<std::size_t>(limit + 1), 0);
    std::vector<int> p(static_cast<std::size_t>(limit / 2 + 1), 0);
    q_odd[0] = 1;
    p[0] = 1;

    for (int k = 1; k <= kmax; ++k) {
        const int k_sq = k * k;
        const int rem = limit - k_sq;
        int idx = k_sq - 2;

        const int vmax_small = std::min(k - 1, rem / 2);
        for (int v = 0; v <= vmax_small; ++v) {
            idx += 2;
            const int pos = k_sq + (v << 1);
            int x = q_odd[static_cast<std::size_t>(pos)] + p[static_cast<std::size_t>(v)];
            if (x >= MOD) x -= MOD;
            q_odd[static_cast<std::size_t>(pos)] = x;
        }

        const int vmax = rem / 2;
        for (int v = k; v <= vmax; ++v) {
            idx += 2;
            int pv = p[static_cast<std::size_t>(v)] + p[static_cast<std::size_t>(v - k)];
            if (pv >= MOD) pv -= MOD;
            p[static_cast<std::size_t>(v)] = pv;
            int x = q_odd[static_cast<std::size_t>(idx)] + pv;
            if (x >= MOD) x -= MOD;
            q_odd[static_cast<std::size_t>(idx)] = x;
        }
    }

    return q_odd;
}

static int solve614(int limit) {
    std::vector<int> q_odd = self_conjugate(limit);
    std::vector<int> s_odd = q_odd;
    for (int i = 1; i <= limit; ++i) {
        int x = s_odd[static_cast<std::size_t>(i)] + s_odd[static_cast<std::size_t>(i - 1)];
        if (x >= MOD) x -= MOD;
        s_odd[static_cast<std::size_t>(i)] = x;
    }

    std::vector<int> q = distinct_parts(limit / 4);

    int res = MOD - 1;
    for (int i = 0; i <= limit; i += 4) {
        res = static_cast<int>((res + 1LL * s_odd[static_cast<std::size_t>(limit - i)] *
                                           q[static_cast<std::size_t>(i / 4)]) %
                               MOD);
    }
    return res;
}

static int brute_cumulative(int limit) {
    std::vector<int> dp(static_cast<std::size_t>(limit + 1), 0);
    dp[0] = 1;
    for (int p = 1; p <= limit; ++p) {
        if ((p & 3) == 2) continue;
        for (int s = limit; s >= p; --s) {
            int x = dp[static_cast<std::size_t>(s)] + dp[static_cast<std::size_t>(s - p)];
            if (x >= MOD) x -= MOD;
            dp[static_cast<std::size_t>(s)] = x;
        }
    }
    int sum = 0;
    for (int s = 1; s <= limit; ++s) {
        sum += dp[static_cast<std::size_t>(s)];
        if (sum >= MOD) sum -= MOD;
    }
    return sum;
}

int main(int argc, char** argv) {
    int n = 10'000'000;
    if (argc >= 2) n = std::stoi(argv[1]);

    assert(solve614(100) == brute_cumulative(100));
    assert(solve614(500) == brute_cumulative(500));
    assert(solve614(1000) == brute_cumulative(1000));

    std::cout << solve614(n) << '\n';
    return 0;
}
