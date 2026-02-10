#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::vector<u32> primes_up_to(u32 n) {
    std::vector<u32> primes;
    if (n < 2) return primes;
    primes.push_back(2);

    const u32 m = n / 2;
    std::vector<std::uint8_t> is_comp(static_cast<std::size_t>(m + 1), 0);

    const u32 r = static_cast<u32>(std::sqrt(static_cast<long double>(n)));
    for (u32 i = 1; 2 * i + 1 <= r; ++i) {
        if (is_comp[i]) continue;
        const u32 p = 2 * i + 1;
        u64 start = (static_cast<u64>(p) * p - 1) / 2;
        for (u64 j = start; j <= m; j += p) is_comp[static_cast<std::size_t>(j)] = 1;
    }

    for (u32 i = 1; i <= m; ++i) {
        if (!is_comp[i]) primes.push_back(2 * i + 1);
    }
    return primes;
}

static u64 vp_fact(u64 n, u32 p) {
    u64 s = 0;
    while (n) {
        n /= p;
        s += n;
    }
    return s;
}

static u32 min_m_for_prime_power(u32 p, u32 exp) {
    // Smallest m such that v_p(m!) >= exp.
    u32 lo = 1;
    u32 hi = p * exp;  // always sufficient since v_p((p*exp)!) >= exp.
    while (lo < hi) {
        const u32 mid = lo + (hi - lo) / 2;
        if (vp_fact(mid, p) >= exp) hi = mid;
        else lo = mid + 1;
    }
    return lo;
}

static u64 solve(u32 n) {
    std::vector<u32> s(static_cast<std::size_t>(n + 1), 0U);
    const std::vector<u32> primes = primes_up_to(n);

    for (u32 p : primes) {
        // exp=1: s(m) >= p for all multiples of p (and s(p)=p).
        for (u32 k = p; k <= n; k += p) {
            if (s[k] < p) s[k] = p;
        }

        // exp>=2: update multiples of p^exp with the required minimal m.
        u64 q = static_cast<u64>(p) * static_cast<u64>(p);
        u32 exp = 2;
        while (q <= n) {
            const u32 qq = static_cast<u32>(q);
            const u32 m = min_m_for_prime_power(p, exp);
            for (u32 k = qq; k <= n; k += qq) {
                if (s[k] < m) s[k] = m;
            }
            if (q > static_cast<u64>(n) / p) break;
            q *= p;
            ++exp;
        }
    }

    u64 sum = 0;
    for (u32 i = 2; i <= n; ++i) sum += static_cast<u64>(s[i]);
    return sum;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    assert(solve(100) == 2012ULL);

    std::cout << solve(100'000'000U) << '\n';
    return 0;
}
