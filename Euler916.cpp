#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 kMod = 1'000'000'007LL;

i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    i64 x = a % kMod;
    i64 p = e;
    while (p > 0) {
        if (p & 1LL) {
            r = static_cast<i64>((__int128)r * x % kMod);
        }
        x = static_cast<i64>((__int128)x * x % kMod);
        p >>= 1LL;
    }
    return r;
}

i64 mod_inv(i64 a) {
    return mod_pow((a % kMod + kMod) % kMod, kMod - 2);
}

i64 P_formula(std::int64_t n) {
    const std::int64_t two_n = 2 * n;

    i64 fact_n = 1;
    i64 fact_2n = 1;
    for (std::int64_t i = 1; i <= two_n; ++i) {
        fact_2n = static_cast<i64>((__int128)fact_2n * (i % kMod) % kMod);
        if (i == n) {
            fact_n = fact_2n;
        }
    }

    const i64 inv_fact_n = mod_inv(fact_n);
    const i64 c2n_n = static_cast<i64>((__int128)fact_2n * inv_fact_n % kMod * inv_fact_n % kMod);

    const i64 inv_n1 = mod_inv(n + 1);
    const i64 inv_n2 = mod_inv(n + 2);

    const i64 f1 = static_cast<i64>((__int128)c2n_n * inv_n1 % kMod);
    const i64 c2n_n1 = static_cast<i64>((__int128)c2n_n * (n % kMod) % kMod * inv_n1 % kMod);
    const i64 f2 = static_cast<i64>((__int128)3 * c2n_n1 % kMod * inv_n2 % kMod);

    i64 ans = static_cast<i64>((__int128)f1 * f1 % kMod);
    ans += static_cast<i64>((__int128)f2 * f2 % kMod);
    if (ans >= kMod) {
        ans -= kMod;
    }
    return ans;
}

int lis_len(const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    std::vector<int> dp(n, 1);
    int best = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (p[j] < p[i]) {
                dp[i] = std::max(dp[i], dp[j] + 1);
            }
        }
        best = std::max(best, dp[i]);
    }
    return best;
}

int lds_len(const std::vector<int>& p) {
    const int n = static_cast<int>(p.size());
    std::vector<int> dp(n, 1);
    int best = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (p[j] > p[i]) {
                dp[i] = std::max(dp[i], dp[j] + 1);
            }
        }
        best = std::max(best, dp[i]);
    }
    return best;
}

std::int64_t brute_count(int n) {
    const int m = 2 * n;
    std::vector<int> p(m);
    for (int i = 0; i < m; ++i) {
        p[i] = i + 1;
    }

    std::int64_t cnt = 0;
    do {
        if (lis_len(p) <= n + 1 && lds_len(p) <= 2) {
            ++cnt;
        }
    } while (std::next_permutation(p.begin(), p.end()));

    return cnt;
}

void validate() {
    for (int n = 1; n <= 4; ++n) {
        assert(P_formula(n) == brute_count(n));
    }

    assert(P_formula(2) == 13);
    assert(P_formula(10) == 45'265'702);
}

}  // namespace

int main() {
    validate();
    std::cout << P_formula(100'000'000) << '\n';
    return 0;
}
