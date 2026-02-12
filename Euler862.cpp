#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u128 S_value(int k) {
    std::array<u128, 13> fact{};
    fact[0] = 1;
    for (int i = 1; i <= k; ++i) fact[i] = fact[i - 1] * static_cast<u128>(i);

    std::array<int, 10> cnt{};
    u128 ans = 0;

    auto dfs = [&](auto&& self, int d, int rem) -> void {
        if (d == 9) {
            cnt[9] = rem;
            u128 perms = fact[k];
            for (int i = 0; i < 10; ++i) perms /= fact[cnt[i]];
            u128 valid = perms * static_cast<u128>(k - cnt[0]) / static_cast<u128>(k);
            ans += valid * (valid - 1) / 2;
            return;
        }

        for (int v = 0; v <= rem; ++v) {
            cnt[d] = v;
            self(self, d + 1, rem - v);
        }
    };

    dfs(dfs, 0, k);
    return ans;
}

static void print_u128(u128 x) {
    if (x == 0) {
        std::cout << 0 << '\n';
        return;
    }
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    for (auto it = s.rbegin(); it != s.rend(); ++it) std::cout << *it;
    std::cout << '\n';
}

int main() {
    assert(S_value(3) == 1701);
    print_u128(S_value(12));
    return 0;
}
