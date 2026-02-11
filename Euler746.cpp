#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kMod = 1'000'000'007U;

u32 add_mod(const u32 a, const u32 b) {
    const u32 s = a + b;
    return (s >= kMod) ? (s - kMod) : s;
}

u32 sub_mod(const u32 a, const u32 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

u32 mul_mod(const u64 a, const u64 b) {
    return static_cast<u32>((a * b) % kMod);
}

u32 mod_pow(u32 base, u64 exp) {
    u32 result = 1U;
    while (exp > 0U) {
        if ((exp & 1U) != 0U) {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exp >>= 1U;
    }
    return result;
}

std::vector<std::vector<u32>> build_cycle_pattern_counts(const int n_max) {
    const int max_len = 2 * n_max;
    const int max_k = n_max;
    std::vector<std::vector<u32>> counts(max_len + 1, std::vector<u32>(max_k + 1, 0U));

    struct Init {
        int a;
        int b;
    };
    std::vector<Init> inits;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            if (a != 0 && b != 0) {
                continue;
            }
            inits.push_back({a, b});
        }
    }

    constexpr int state_count = 9;
    for (const Init init : inits) {
        std::vector<std::vector<u32>> dp(state_count, std::vector<u32>(max_k + 1, 0U));
        std::vector<std::vector<u32>> next_dp(state_count, std::vector<u32>(max_k + 1, 0U));

        const int k0 = (init.a != 0) + (init.b != 0);
        const int state0 = 3 * init.a + init.b;
        dp[state0][k0] = 1U;

        for (int len = 2; len <= max_len; ++len) {
            for (int state = 0; state < state_count; ++state) {
                const int prev2 = state / 3;
                const int prev1 = state % 3;
                if (prev1 != 0 && init.a != 0) {
                    continue;
                }
                if (prev2 == 2 && init.a == 1) {
                    continue;
                }
                if (prev1 == 2 && init.b == 1) {
                    continue;
                }

                const int k_cap = std::min(max_k, len / 2);
                for (int k = 0; k <= k_cap; ++k) {
                    const u32 ways = dp[state][k];
                    if (ways == 0U) {
                        continue;
                    }
                    counts[len][k] = add_mod(counts[len][k], ways);
                }
            }

            if (len == max_len) {
                break;
            }

            for (int state = 0; state < state_count; ++state) {
                std::fill(next_dp[state].begin(), next_dp[state].end(), 0U);
            }

            const int k_cap = std::min(max_k, (len + 1) / 2);
            for (int state = 0; state < state_count; ++state) {
                const int prev2 = state / 3;
                const int prev1 = state % 3;
                for (int k = 0; k <= k_cap; ++k) {
                    const u32 ways = dp[state][k];
                    if (ways == 0U) {
                        continue;
                    }
                    for (int cur = 0; cur < 3; ++cur) {
                        if (prev1 != 0 && cur != 0) {
                            continue;
                        }
                        if (prev2 == 2 && cur == 1) {
                            continue;
                        }
                        const int nk = k + (cur != 0);
                        if (nk > max_k) {
                            continue;
                        }
                        const int next_state = 3 * prev1 + cur;
                        next_dp[next_state][nk] = add_mod(next_dp[next_state][nk], ways);
                    }
                }
            }

            dp.swap(next_dp);
        }
    }

    return counts;
}

std::vector<u32> compute_m_values(const int n_max) {
    const auto counts = build_cycle_pattern_counts(n_max);
    const int max_fact = 2 * n_max;

    std::vector<u32> fact(max_fact + 1, 1U);
    for (int i = 1; i <= max_fact; ++i) {
        fact[i] = mul_mod(fact[i - 1], static_cast<u32>(i));
    }

    std::vector<u32> inv_fact(max_fact + 1, 1U);
    inv_fact[max_fact] = mod_pow(fact[max_fact], kMod - 2U);
    for (int i = max_fact; i >= 1; --i) {
        inv_fact[i - 1] = mul_mod(inv_fact[i], static_cast<u32>(i));
    }

    std::vector<u32> pow4(n_max + 1, 1U);
    for (int i = 1; i <= n_max; ++i) {
        pow4[i] = mul_mod(pow4[i - 1], 4U);
    }

    const auto falling = [&](const int n, const int k) -> u32 {
        return mul_mod(fact[n], inv_fact[n - k]);
    };

    std::vector<u32> m_values(n_max + 1, 0U);
    for (int n = 2; n <= n_max; ++n) {
        u32 total = 0U;
        for (int k = 0; k <= n; ++k) {
            const u32 pattern_ways = counts[2 * n][k];
            if (pattern_ways == 0U) {
                continue;
            }
            const int rem = 2 * n - 2 * k;
            u32 term = pattern_ways;
            term = mul_mod(term, falling(n, k));
            term = mul_mod(term, pow4[k]);
            term = mul_mod(term, mul_mod(fact[rem], fact[rem]));
            total = (k % 2 == 0) ? add_mod(total, term) : sub_mod(total, term);
        }
        m_values[n] = mul_mod(2U, total);
    }

    return m_values;
}

u32 compute_s_value(const int n) {
    const auto m_values = compute_m_values(n);
    u32 total = 0U;
    for (int k = 2; k <= n; ++k) {
        total = add_mod(total, m_values[k]);
    }

    assert(m_values[1] == 0U);
    assert(m_values[2] == 896U);
    assert(m_values[3] == 890'880U);
    assert(m_values[10] == 170'717'180U);

    u32 s10 = 0U;
    for (int k = 2; k <= 10; ++k) {
        s10 = add_mod(s10, m_values[k]);
    }
    assert(s10 == 399'291'975U);

    return total;
}

}  // namespace

int main() {
    std::cout << compute_s_value(2021) << '\n';
    return 0;
}
