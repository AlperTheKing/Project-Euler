#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 clmul(u64 x, u64 y) {
    u64 r = 0;
    while (y > 0) {
        if (y & 1ULL) r ^= x;
        x <<= 1ULL;
        y >>= 1ULL;
    }
    return r;
}

static bool is_solution(u64 a, u64 b) {
    u64 lhs = clmul(a, a) ^ clmul(clmul(2, a), b) ^ clmul(b, b);
    return lhs == 5;
}

static std::vector<std::pair<u64, u64>> generated_solutions(u64 n) {
    std::vector<std::pair<u64, u64>> sol;
    u64 a = 0;
    u64 b = 3;
    while (b <= n) {
        sol.push_back({a, b});
        u64 next = static_cast<u64>((static_cast<u128>(b) << 1) ^ static_cast<u128>(a));
        a = b;
        b = next;
    }
    return sol;
}

static std::vector<std::pair<u64, u64>> brute_solutions(u64 n) {
    std::vector<std::pair<u64, u64>> sol;
    for (u64 b = 0; b <= n; ++b) {
        for (u64 a = 0; a <= b; ++a) {
            if (is_solution(a, b)) sol.push_back({a, b});
        }
    }
    return sol;
}

static u64 x_fast(u64 n) {
    u64 x = 0;
    u64 a = 0;
    u64 b = 3;
    while (b <= n) {
        x ^= b;
        u64 next = static_cast<u64>((static_cast<u128>(b) << 1) ^ static_cast<u128>(a));
        a = b;
        b = next;
    }
    return x;
}

static u64 x_bruteforce(u64 n) {
    u64 x = 0;
    for (u64 b = 0; b <= n; ++b) {
        for (u64 a = 0; a <= b; ++a) {
            if (is_solution(a, b)) x ^= b;
        }
    }
    return x;
}

int main() {
    assert(x_fast(10) == 5ULL);

    {
        auto g = generated_solutions(2000);
        auto b = brute_solutions(2000);
        assert(g == b);
        assert(x_fast(2000) == x_bruteforce(2000));
    }

    std::cout << x_fast(1'000'000'000'000'000'000ULL) << '\n';
    return 0;
}
