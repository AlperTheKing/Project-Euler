#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 MOD = 1'000'000'007LL;

i64 norm(i64 v) {
    v %= MOD;
    if (v < 0) {
        v += MOD;
    }
    return v;
}

std::vector<i64> brute_values(int n) {
    std::vector<std::vector<i64>> cnt(
        static_cast<std::size_t>(n / 2 + 2),
        std::vector<i64>(static_cast<std::size_t>(n + 1), 0LL));

    for (int a = 0; 2 * a <= n; ++a) {
        for (int b = 0; 2 * a + 3 * b <= n; ++b) {
            for (int c = 0; 2 * a + 3 * b + 5 * c <= n; ++c) {
                const int s = 2 * a + 3 * b + 5 * c;
                const int k = a + b + c;
                cnt[static_cast<std::size_t>(k)][static_cast<std::size_t>(s)] += 1LL;
            }
        }
    }

    std::vector<i64> out(static_cast<std::size_t>(n + 1), 0LL);
    for (int t = 0; t <= n; ++t) {
        i64 val = 0LL;
        for (std::size_t k = 0; k < cnt.size(); ++k) {
            for (int s = 0; s <= t; ++s) {
                val += cnt[k][static_cast<std::size_t>(s)] *
                       cnt[k][static_cast<std::size_t>(t - s)];
            }
        }
        out[static_cast<std::size_t>(t)] = val;
    }

    return out;
}

i64 fast_f(int n) {
    std::vector<i64> a(static_cast<std::size_t>(n + 1), 0LL);
    std::vector<i64> num(13, 0LL);
    num[0] = 1;
    num[1] = -1;
    num[4] = 1;
    num[5] = 1;
    num[6] = -2;
    num[7] = 1;
    num[8] = 1;
    num[11] = -1;
    num[12] = 1;

    for (int i = 0; i <= n; ++i) {
        i64 v = (i <= 12) ? num[static_cast<std::size_t>(i)] : 0LL;

        if (i >= 1) v += a[static_cast<std::size_t>(i - 1)];
        if (i >= 6) v += a[static_cast<std::size_t>(i - 6)];
        if (i >= 9) v -= a[static_cast<std::size_t>(i - 9)];
        if (i >= 10) v += a[static_cast<std::size_t>(i - 10)];
        if (i >= 11) v -= a[static_cast<std::size_t>(i - 11)];
        if (i >= 13) v -= a[static_cast<std::size_t>(i - 13)];
        if (i >= 19) v += a[static_cast<std::size_t>(i - 19)];
        if (i >= 21) v += a[static_cast<std::size_t>(i - 21)];
        if (i >= 22) v -= a[static_cast<std::size_t>(i - 22)];
        if (i >= 23) v += a[static_cast<std::size_t>(i - 23)];
        if (i >= 26) v -= a[static_cast<std::size_t>(i - 26)];
        if (i >= 31) v -= a[static_cast<std::size_t>(i - 31)];
        if (i >= 32) v += a[static_cast<std::size_t>(i - 32)];

        a[static_cast<std::size_t>(i)] = norm(v);
    }

    return a[static_cast<std::size_t>(n)];
}

}  // namespace

int main() {
    const std::vector<i64> brute = brute_values(120);
    assert(brute[10] == 4LL);
    assert(brute[100] == 3'629LL);

    assert(fast_f(10) == 4LL);
    assert(fast_f(100) == 3'629LL);

    for (int i = 0; i <= 120; ++i) {
        assert(fast_f(i) == brute[static_cast<std::size_t>(i)] % MOD);
    }

    std::cout << fast_f(10'000'000) << "\n";
    return 0;
}
