#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'020'202'009ULL;

u64 f_bruteforce(const int n) {
    std::vector<std::vector<u64>> comb(static_cast<std::size_t>(n + 1),
                                       std::vector<u64>(static_cast<std::size_t>(n + 1), 0ULL));
    for (int i = 0; i <= n; ++i) {
        comb[static_cast<std::size_t>(i)][0] = 1ULL;
        comb[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1ULL;
        for (int j = 1; j < i; ++j) {
            comb[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                comb[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)] +
                comb[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)];
        }
    }

    std::vector<u64> dp(static_cast<std::size_t>(n + 2), 0ULL);
    dp[1] = 1ULL;
    for (int k = 2; k <= n; ++k) {
        std::vector<u64> next(static_cast<std::size_t>(n + 2), 0ULL);
        for (int r_old = 1; r_old <= n; ++r_old) {
            const u64 ways = dp[static_cast<std::size_t>(r_old)];
            if (ways == 0ULL) {
                continue;
            }
            for (int t = 0; t <= r_old; t += 2) {
                const int r_new = r_old - t + 1;
                next[static_cast<std::size_t>(r_new)] +=
                    ways * comb[static_cast<std::size_t>(r_old)][static_cast<std::size_t>(t)];
            }
        }
        dp.swap(next);
    }

    u64 total = 0ULL;
    for (const u64 v : dp) {
        total += v;
    }
    return total;
}

u64 solve(const int n) {
    if (n == 0) {
        return 1ULL;
    }
    std::vector<u32> prev(static_cast<std::size_t>(n + 1), 0U);
    std::vector<u32> curr(static_cast<std::size_t>(n + 1), 0U);
    prev[0] = 1U;

    for (int row = 1; row <= n; ++row) {
        curr[0] = 0U;
        for (int k = 1; k <= row; ++k) {
            u32 v = curr[static_cast<std::size_t>(k - 1)] +
                    prev[static_cast<std::size_t>(row - k)];
            if (v >= kMod) {
                v -= static_cast<u32>(kMod);
            }
            curr[static_cast<std::size_t>(k)] = v;
        }
        std::swap(prev, curr);
    }
    return prev[static_cast<std::size_t>(n)];
}

}  // namespace

int main() {
    assert(f_bruteforce(4) == 5ULL);
    assert(f_bruteforce(8) == 1385ULL);
    assert(solve(8) == 1385ULL);
    assert(solve(12) == f_bruteforce(12) % kMod);

    std::cout << solve(24'680) << '\n';
    return 0;
}
