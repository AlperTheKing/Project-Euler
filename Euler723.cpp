#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::string to_string_u128(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v > 0) {
        const int d = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + d));
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

template <typename F>
static u128 product_over_exps(const std::array<int, 8>& exps, F term) {
    u128 out = 1;
    for (const int e : exps) out *= static_cast<u128>(term(e));
    return out;
}

static u128 F_value(const std::array<int, 8>& exps) {
    const u128 t1 = 7 * product_over_exps(exps, [](int e) -> u64 {
        return static_cast<u64>(e + 1);
    });
    const u128 t2 = 14 * product_over_exps(exps, [](int e) -> u64 {
        const u64 x = static_cast<u64>(e + 1);
        return x * x;
    });
    const u128 t3 = 4 * product_over_exps(exps, [](int e) -> u64 {
        const u64 x = static_cast<u64>(e + 1);
        return (x * x + 1) / 2;
    });
    const u128 t4 = 8 * product_over_exps(exps, [](int e) -> u64 {
        const u64 x = static_cast<u64>(e + 1);
        return x * x * x;
    });
    const u128 t5 = 4 * product_over_exps(exps, [](int e) -> u64 {
        const u64 x = static_cast<u64>(e + 1);
        const u64 y = static_cast<u64>(2 * e * e + 4 * e + 3);
        return (x * y) / 3;
    });
    return t1 - t2 - t3 + t4 + t5;
}

static u128 S_value(const std::array<int, 8>& exps) {
    const u128 t1 = 7 * product_over_exps(exps, [](int e) -> u64 {
        return static_cast<u64>((e + 1) * (e + 2) / 2);
    });
    const u128 t2 = 14 * product_over_exps(exps, [](int e) -> u64 {
        return static_cast<u64>((e + 1) * (e + 2) * (2 * e + 3) / 6);
    });
    const u128 t3 = 4 * product_over_exps(exps, [](int e) -> u64 {
        const u64 a = static_cast<u64>((e + 1) * (e + 2) * (2 * e + 3) / 6);
        const u64 b = static_cast<u64>(e / 2 + 1);
        return (a + b) / 2;
    });
    const u128 t4 = 8 * product_over_exps(exps, [](int e) -> u64 {
        const u64 x = static_cast<u64>((e + 1) * (e + 2) / 2);
        return x * x;
    });
    const u128 t5 = 4 * product_over_exps(exps, [](int e) -> u64 {
        return static_cast<u64>((e + 1) * (e + 2) * (e * e + 3 * e + 3) / 6);
    });
    return t1 - t2 - t3 + t4 + t5;
}

int main() {
    {
        std::array<int, 8> exps{};
        assert(F_value(exps) == 1);
    }
    {
        std::array<int, 8> exps{};
        exps[0] = 1;
        assert(F_value(exps) == 38);
    }
    {
        std::array<int, 8> exps{};
        exps[0] = 2;
        assert(F_value(exps) == 167);
    }
    {
        std::array<int, 8> exps{};
        exps[0] = 2;
        exps[1] = 1;
        assert(S_value(exps) == 2370);
    }
    {
        std::array<int, 8> exps{};
        exps[0] = 1;
        exps[1] = 1;
        exps[2] = 1;
        assert(S_value(exps) == 5535);
    }

    const std::array<int, 8> target_exps = {6, 3, 2, 1, 1, 1, 1, 1};
    std::cout << to_string_u128(S_value(target_exps)) << '\n';
    return 0;
}
