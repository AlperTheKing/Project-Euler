#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr u64 kMod = 1'020'202'009ULL;

u64 mod_pow(u64 a, u64 e) {
    u64 r = 1ULL;
    while (e > 0ULL) {
        if (e & 1ULL) {
            r = static_cast<u64>((static_cast<unsigned __int128>(r) * a) % kMod);
        }
        a = static_cast<u64>((static_cast<unsigned __int128>(a) * a) % kMod);
        e >>= 1ULL;
    }
    return r;
}

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
    std::vector<u64> fact(static_cast<std::size_t>(n + 1), 1ULL);
    std::vector<u64> inv_fact(static_cast<std::size_t>(n + 1), 1ULL);

    for (int i = 1; i <= n; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<u64>((static_cast<unsigned __int128>(
                                  fact[static_cast<std::size_t>(i - 1)]) *
                              static_cast<u64>(i)) %
                             kMod);
    }
    inv_fact[static_cast<std::size_t>(n)] = mod_pow(fact[static_cast<std::size_t>(n)], kMod - 2ULL);
    for (int i = n; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<u64>((static_cast<unsigned __int128>(
                                  inv_fact[static_cast<std::size_t>(i)]) *
                              static_cast<u64>(i)) %
                             kMod);
    }

    auto binom = [&](const int a, const int b) -> u64 {
        if (b < 0 || b > a) {
            return 0ULL;
        }
        u64 v = fact[static_cast<std::size_t>(a)];
        v = static_cast<u64>((static_cast<unsigned __int128>(v) *
                              inv_fact[static_cast<std::size_t>(b)]) %
                             kMod);
        v = static_cast<u64>((static_cast<unsigned __int128>(v) *
                              inv_fact[static_cast<std::size_t>(a - b)]) %
                             kMod);
        return v;
    };

    const u64 inv2 = (kMod + 1ULL) / 2ULL;
    std::vector<u64> E(static_cast<std::size_t>(n + 1), 0ULL);
    E[0] = 1ULL;
    if (n >= 1) {
        E[1] = 1ULL;
    }

    for (int m = 1; m < n; ++m) {
        u64 sum = 0ULL;
        for (int k = 0; k <= m; ++k) {
            u64 term = binom(m, k);
            term = static_cast<u64>((static_cast<unsigned __int128>(term) *
                                     E[static_cast<std::size_t>(k)]) %
                                    kMod);
            term = static_cast<u64>((static_cast<unsigned __int128>(term) *
                                     E[static_cast<std::size_t>(m - k)]) %
                                    kMod);
            sum += term;
            if (sum >= kMod) {
                sum -= kMod;
            }
        }
        E[static_cast<std::size_t>(m + 1)] =
            static_cast<u64>((static_cast<unsigned __int128>(sum) * inv2) % kMod);
    }

    return E[static_cast<std::size_t>(n)];
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
