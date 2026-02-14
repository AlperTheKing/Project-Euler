#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 cube_u64(u64 k) {
    return static_cast<u64>(static_cast<u128>(k) * k * k);
}

static u64 floor_cuberoot(u64 x) {
    long double y = std::cbrt(static_cast<long double>(x));
    u64 k = static_cast<u64>(y);
    while (cube_u64(k + 1) <= x) ++k;
    while (cube_u64(k) > x) --k;
    return k;
}

static u64 D_bruteforce(u64 n) {
    u64 steps = 0;
    while (n > 0) {
        u64 k = floor_cuberoot(n);
        n -= cube_u64(k);
        ++steps;
    }
    return steps;
}

static u128 S_bruteforce(u64 N) {
    u128 s = 0;
    for (u64 n = 1; n < N; ++n) s += D_bruteforce(n);
    return s;
}

static u128 S_fast(u64 N) {
    if (N <= 1) return 0;

    u64 K = floor_cuberoot(N - 1);
    std::vector<u128> A(static_cast<std::size_t>(K + 1), 0);

    std::unordered_map<u64, u128> memo;
    memo.reserve(static_cast<std::size_t>(K * 3));
    memo[0] = 0;
    memo[1] = 0;
    if (K >= 1) A[1] = 0;

    std::function<u128(u64, u64)> solve = [&](u64 x, u64 ready_k) -> u128 {
        if (x <= 1) return 0;
        auto it = memo.find(x);
        if (it != memo.end()) return it->second;

        u64 k = floor_cuberoot(x - 1);
        assert(k <= ready_k);
        u64 c = cube_u64(k);
        u64 t = x - c;

        u128 ans = A[static_cast<std::size_t>(k)] + static_cast<u128>(t) + solve(t, ready_k);
        memo[x] = ans;
        return ans;
    };

    for (u64 k = 1; k < K; ++k) {
        u64 delta = 3 * k * k + 3 * k + 1;
        u128 sdelta = solve(delta, k);
        A[static_cast<std::size_t>(k + 1)] = A[static_cast<std::size_t>(k)] + static_cast<u128>(delta) + sdelta;
        memo[cube_u64(k + 1)] = A[static_cast<std::size_t>(k + 1)];
    }

    return solve(N, K);
}

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
    assert(S_fast(100) == 512);
    for (u64 n : {2ULL, 10ULL, 100ULL, 500ULL, 1000ULL}) {
        assert(S_fast(n) == S_bruteforce(n));
    }

    std::cout << to_string_u128(S_fast(100'000'000'000'000'000ULL)) << '\n';
    return 0;
}
