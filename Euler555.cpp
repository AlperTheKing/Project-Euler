#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

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

static u128 S(u64 p, u64 m) {
    u128 total = 0;
    for (u64 g = 1; g <= p; ++g) {
        const u64 T = p / g;
        if (T < 2) continue;  // need k = g*t with t>=2 (i.e., g<k)

        const u64 cnt = T - 1;                              // t = 2..T
        const u128 sum_t = (u128)T * (T + 1) / 2 - 1;       // sum_{t=2..T} t

        const u128 A1 = (u128)cnt * g * (m + g + 1);        // cnt * g * (m+g+1)
        const u128 A2 = (u128)cnt * g * (g - 1) / 2;        // cnt * g*(g-1)/2
        const u128 B = (u128)g * g * sum_t;                 // g^2 * sum_t

        total += (A1 + A2 - B);
    }
    return total;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    assert(to_string_u128(S(10, 10)) == "225");
    assert(to_string_u128(S(1000, 1000)) == "208724467");

    std::cout << to_string_u128(S(1'000'000ULL, 1'000'000ULL)) << '\n';
    return 0;
}

