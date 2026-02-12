#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::pair<u64, u64> fib_pair(u64 n, u64 mod) {
    if (n == 0) return {0, 1 % mod};
    auto [a, b] = fib_pair(n >> 1, mod);
    u64 c = static_cast<u64>((static_cast<u128>(a) * ((2 * static_cast<u128>(b) + mod - a) % mod)) % mod);
    u64 d = static_cast<u64>((static_cast<u128>(a) * a + static_cast<u128>(b) * b) % mod);
    if ((n & 1ULL) == 0) return {c, d};
    return {d, (c + d) % mod};
}

static std::vector<u64> divisors_of(u64 n) {
    std::vector<u64> divs;
    for (u64 d = 1; d * d <= n; ++d) {
        if (n % d == 0) {
            divs.push_back(d);
            if (d * d != n) divs.push_back(n / d);
        }
    }
    std::sort(divs.begin(), divs.end());
    return divs;
}

static bool has_period(u64 mod, u64 period, const std::vector<u64>& proper_divs) {
    auto p = fib_pair(period, mod);
    if (p.first != 0 || p.second != 1 % mod) return false;
    for (u64 d : proper_divs) {
        auto q = fib_pair(d, mod);
        if (q.first == 0 && q.second == 1 % mod) return false;
    }
    return true;
}

static u64 sample_check() {
    std::vector<u64> proper;
    for (u64 d : divisors_of(18)) {
        if (d < 18) proper.push_back(d);
    }
    u64 sum = 0;
    for (u64 n = 2; n < 50; ++n) {
        if (has_period(n, 18, proper)) sum += n;
    }
    return sum;
}

int main() {
    if (sample_check() != 57) return 1;

    const u64 limit = 1'000'000'000ULL;
    const u64 period = 120;
    std::vector<u64> proper_divs;
    for (u64 d : divisors_of(period)) {
        if (d < period) proper_divs.push_back(d);
    }

    const std::vector<std::pair<u64, int>> fac = {
        {2, 5}, {3, 2}, {5, 1}, {7, 1}, {11, 1}, {23, 1}, {31, 1},
        {41, 1}, {61, 1}, {241, 1}, {2161, 1}, {2521, 1}, {20641, 1}
    };

    u64 ans = 0;
    std::vector<u64> divisors{1};
    for (auto [p, e] : fac) {
        std::vector<u64> next;
        next.reserve(divisors.size() * (e + 1));
        u64 pe = 1;
        for (int i = 0; i <= e; ++i) {
            for (u64 d : divisors) {
                u128 v = static_cast<u128>(d) * pe;
                if (v < limit) next.push_back(static_cast<u64>(v));
            }
            pe *= p;
        }
        divisors.swap(next);
    }

    for (u64 n : divisors) {
        if (n < 2) continue;
        if (has_period(n, period, proper_divs)) ans += n;
    }

    std::cout << ans << '\n';
    return 0;
}
