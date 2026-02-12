#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct BoundedSearch {
    u64 limit;
    std::vector<std::unordered_map<std::uint32_t, u64>> memo;

    explicit BoundedSearch(u64 lim) : limit(lim), memo(128) {
        for (auto& mp : memo) mp.reserve(4096);
    }

    u64 F(int m, long long d) {
        if (d < -1) return 0;
        if (d > static_cast<long long>(limit)) d = static_cast<long long>(limit);
        if (m == 0) return 1;
        if (d == -1) return 1;
        if (d >= m - 1) {
            if (m >= 63) return limit;
            u64 full = 1ULL << m;
            return std::min(full, limit);
        }

        std::uint32_t key = static_cast<std::uint32_t>(d + 1);
        auto& mp = memo[m];
        auto it = mp.find(key);
        if (it != mp.end()) return it->second;

        u64 y = F(m - 1, d - 1);
        long long d2 = static_cast<long long>(y) + d - 1;
        u64 z = F(m - 1, d2);

        u64 ans = y + z;
        if (ans > limit) ans = limit;
        mp[key] = ans;
        return ans;
    }

    int Q(u64 N, int d) {
        if (d == 0) return static_cast<int>(N - 1);
        int m = 0;
        while (F(m, d) < N) ++m;
        return m;
    }
};

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    u64 L = 1;
    for (int i = 0; i < 10; ++i) L *= 7ULL;

    BoundedSearch bs(L);

    assert(bs.Q(7, 1) == 3);
    assert(bs.Q(777, 2) == 10);
    for (u64 n = 1; n <= 20; ++n) {
        assert(bs.Q(n, 0) == static_cast<int>(n - 1));
    }

    u128 ans = static_cast<u128>(L) * static_cast<u128>(L - 1) / 2;
    for (int d = 1; d <= 7; ++d) {
        u64 prev = 0;
        int m = 0;
        while (prev < L) {
            u64 cur = bs.F(m, d);
            u64 upto = std::min(cur, L);
            ans += static_cast<u128>(m) * static_cast<u128>(upto - prev);
            prev = upto;
            ++m;
        }
    }

    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
