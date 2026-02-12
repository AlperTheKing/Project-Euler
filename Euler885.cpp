#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 S_mod(int n, u64 mod) {
    std::vector<u64> fact(static_cast<std::size_t>(n + 1), 1);
    for (int i = 1; i <= n; ++i) fact[static_cast<std::size_t>(i)] = fact[static_cast<std::size_t>(i - 1)] * static_cast<u64>(i);

    std::vector<u64> pow10_exact(static_cast<std::size_t>(n + 1), 1);
    std::vector<u64> rep_exact(static_cast<std::size_t>(n + 1), 0);
    std::vector<u64> pow10_mod(static_cast<std::size_t>(n + 1), 1 % mod);
    std::vector<u64> rep_mod(static_cast<std::size_t>(n + 1), 0);

    for (int i = 1; i <= n; ++i) {
        pow10_exact[static_cast<std::size_t>(i)] = pow10_exact[static_cast<std::size_t>(i - 1)] * 10ULL;
        pow10_mod[static_cast<std::size_t>(i)] = (pow10_mod[static_cast<std::size_t>(i - 1)] * 10ULL) % mod;
    }
    for (int i = 0; i <= n; ++i) {
        rep_exact[static_cast<std::size_t>(i)] = (pow10_exact[static_cast<std::size_t>(i)] - 1ULL) / 9ULL;
        if (i == 0) rep_mod[static_cast<std::size_t>(i)] = 0;
        else rep_mod[static_cast<std::size_t>(i)] = (rep_mod[static_cast<std::size_t>(i - 1)] * 10ULL + 1ULL) % mod;
    }

    u64 ans = 0;

    std::function<void(int, int, int, u64, u64)> dfs = [&](int digit, int rem, int s, u64 denom, u64 val_mod) {
        if (digit == 10) {
            u64 ways = fact[static_cast<std::size_t>(n)] / (fact[static_cast<std::size_t>(n - s)] * denom);
            ans = (ans + (ways % mod) * val_mod) % mod;
            return;
        }

        for (int c = 0; c <= rem; ++c) {
            u64 next_denom = denom * fact[static_cast<std::size_t>(c)];
            u64 next_mod = (val_mod * pow10_mod[static_cast<std::size_t>(c)] + (static_cast<u64>(digit) * rep_mod[static_cast<std::size_t>(c)]) % mod) % mod;
            dfs(digit + 1, rem - c, s + c, next_denom, next_mod);
        }
    };

    dfs(1, n, 0, 1ULL, 0ULL);
    return ans;
}

static u128 S_exact_small(int n) {
    std::vector<u64> fact(static_cast<std::size_t>(n + 1), 1);
    for (int i = 1; i <= n; ++i) fact[static_cast<std::size_t>(i)] = fact[static_cast<std::size_t>(i - 1)] * static_cast<u64>(i);

    std::vector<u64> pow10(static_cast<std::size_t>(n + 1), 1);
    std::vector<u64> rep(static_cast<std::size_t>(n + 1), 0);
    for (int i = 1; i <= n; ++i) pow10[static_cast<std::size_t>(i)] = pow10[static_cast<std::size_t>(i - 1)] * 10ULL;
    for (int i = 0; i <= n; ++i) rep[static_cast<std::size_t>(i)] = (pow10[static_cast<std::size_t>(i)] - 1ULL) / 9ULL;

    u128 ans = 0;

    std::function<void(int, int, int, u64, u64)> dfs = [&](int digit, int rem, int s, u64 denom, u64 val) {
        if (digit == 10) {
            u64 ways = fact[static_cast<std::size_t>(n)] / (fact[static_cast<std::size_t>(n - s)] * denom);
            ans += static_cast<u128>(ways) * static_cast<u128>(val);
            return;
        }

        for (int c = 0; c <= rem; ++c) {
            u64 next_denom = denom * fact[static_cast<std::size_t>(c)];
            u64 next_val = val * pow10[static_cast<std::size_t>(c)] + static_cast<u64>(digit) * rep[static_cast<std::size_t>(c)];
            dfs(digit + 1, rem - c, s + c, next_denom, next_val);
        }
    };

    dfs(1, n, 0, 1ULL, 0ULL);
    return ans;
}

int main() {
    constexpr u64 MOD = 1'123'455'689ULL;

    assert(S_exact_small(1) == 45ULL);
    assert(S_exact_small(5) == 1'543'545'675ULL);

    std::cout << S_mod(18, MOD) << '\n';
    return 0;
}
