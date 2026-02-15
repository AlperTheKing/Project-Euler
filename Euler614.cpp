#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using i64 = long long;

static constexpr int MOD = 1'000'000'007;

static std::vector<int> build_a015128(int limit) {
    std::vector<int> a(static_cast<std::size_t>(limit + 1), 0);
    a[0] = 1;

    std::vector<int> squares;
    std::vector<int> signs;
    squares.reserve(static_cast<std::size_t>(std::sqrt(static_cast<double>(limit))) + 2);
    signs.reserve(squares.capacity());
    for (int m = 1; 1LL * m * m <= limit; ++m) {
        squares.push_back(m * m);
        signs.push_back((m & 1) ? 1 : -1);
    }

    for (int n = 1; n <= limit; ++n) {
        i64 val = 0;
        int cnt = 0;
        for (std::size_t i = 0; i < squares.size(); ++i) {
            const int sq = squares[i];
            if (sq > n) break;
            val += (signs[i] > 0) ? a[static_cast<std::size_t>(n - sq)]
                                  : -a[static_cast<std::size_t>(n - sq)];
            if ((++cnt & 63) == 0) val %= MOD;
        }
        val %= MOD;
        if (val < 0) val += MOD;
        a[static_cast<std::size_t>(n)] = static_cast<int>((2LL * val) % MOD);
    }

    return a;
}

static int solve614(int limit) {
    const int max_k = limit / 4;
    std::vector<int> a = build_a015128(max_k);

    std::vector<int> mark(static_cast<std::size_t>(limit + 1), 0);
    mark[0] = 1;
    for (int n = 1;; ++n) {
        const i64 h1 = 1LL * n * (2LL * n - 1);
        if (h1 > limit) break;
        mark[static_cast<std::size_t>(h1)] = 1;
        const i64 h2 = 1LL * n * (2LL * n + 1);
        if (h2 <= limit) mark[static_cast<std::size_t>(h2)] = 1;
    }

    std::vector<int> pref(static_cast<std::size_t>(limit + 1), 0);
    int run = 0;
    for (int i = 0; i <= limit; ++i) {
        run += mark[static_cast<std::size_t>(i)];
        pref[static_cast<std::size_t>(i)] = run;
    }

    i64 ans = 0;
    for (int k = 0; k <= max_k; ++k) {
        ans += 1LL * a[static_cast<std::size_t>(k)] *
               pref[static_cast<std::size_t>(limit - 4 * k)];
        ans %= MOD;
    }
    ans = (ans - 1) % MOD;
    if (ans < 0) ans += MOD;
    return static_cast<int>(ans);
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
    for (int i = 1; i < argc; ++i) {
        n = std::stoi(std::string(argv[i]));
    }

    assert(solve614(100) == brute_cumulative(100));
    assert(solve614(500) == brute_cumulative(500));
    assert(solve614(1000) == brute_cumulative(1000));

    std::cout << solve614(n) << '\n';
    return 0;
}
