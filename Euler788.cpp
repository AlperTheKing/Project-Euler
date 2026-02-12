#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = long long;

static constexpr i64 MOD = 1'000'000'007LL;

static i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    while (e > 0) {
        if (e & 1LL) {
            r = (r * a) % MOD;
        }
        a = (a * a) % MOD;
        e >>= 1LL;
    }
    return r;
}

static i64 comb(int n, int k, const std::vector<i64>& fact, const std::vector<i64>& inv_fact) {
    if (k < 0 || k > n) {
        return 0;
    }
    return (((fact[n] * inv_fact[k]) % MOD) * inv_fact[n - k]) % MOD;
}

static i64 solve_fast(int n) {
    std::vector<i64> fact(n + 1, 1);
    std::vector<i64> inv_fact(n + 1, 1);
    std::vector<i64> pow9(n + 2, 1);

    for (int i = 1; i <= n; ++i) {
        fact[i] = (fact[i - 1] * i) % MOD;
    }
    inv_fact[n] = mod_pow(fact[n], MOD - 2);
    for (int i = n; i >= 1; --i) {
        inv_fact[i - 1] = (inv_fact[i] * i) % MOD;
    }
    for (int i = 1; i <= n + 1; ++i) {
        pow9[i] = (pow9[i - 1] * 9) % MOD;
    }

    i64 ans = 0;
    for (int len = 1; len <= n; ++len) {
        for (int k = len / 2 + 1; k <= len; ++k) {
            const i64 c = comb(len, k, fact, inv_fact);
            const i64 ways = (pow9[len - k + 1] * c) % MOD;
            ans += ways;
            if (ans >= MOD) {
                ans -= MOD;
            }
        }
    }

    return ans;
}

static bool is_dominating(i64 x) {
    std::array<int, 10> cnt{};
    int len = 0;
    while (x > 0) {
        ++cnt[static_cast<int>(x % 10)];
        x /= 10;
        ++len;
    }
    const int best = *std::max_element(cnt.begin(), cnt.end());
    return 2 * best > len;
}

static i64 solve_bruteforce(int n) {
    i64 lim = 1;
    for (int i = 0; i < n; ++i) {
        lim *= 10;
    }

    i64 ans = 0;
    for (i64 x = 1; x < lim; ++x) {
        if (is_dominating(x)) {
            ++ans;
        }
    }
    return ans % MOD;
}

int main() {
    assert(solve_fast(4) == 603);
    assert(solve_fast(10) == 21'893'256);
    assert(solve_fast(6) == solve_bruteforce(6));

    std::cout << solve_fast(2022) << '\n';
    return 0;
}
