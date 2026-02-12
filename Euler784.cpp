#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::vector<int> build_spf(int n) {
    std::vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) {
        spf[i] = i;
    }
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (spf[i] != i) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            if (spf[j] == j) {
                spf[j] = i;
            }
        }
    }
    return spf;
}

static std::vector<std::pair<int, int>> factorize_small(int x, const std::vector<int>& spf) {
    std::vector<std::pair<int, int>> out;
    while (x > 1) {
        const int p = spf[x];
        int c = 0;
        while (x % p == 0) {
            x /= p;
            ++c;
        }
        out.push_back({p, c});
    }
    return out;
}

static std::vector<std::pair<int, int>> merge_factorizations(
    const std::vector<std::pair<int, int>>& a,
    const std::vector<std::pair<int, int>>& b
) {
    std::vector<std::pair<int, int>> f = a;
    for (const auto& [p, c] : b) {
        bool found = false;
        for (auto& [q, d] : f) {
            if (q == p) {
                d += c;
                found = true;
                break;
            }
        }
        if (!found) {
            f.push_back({p, c});
        }
    }
    return f;
}

static std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u128 solve_fast(int n) {
    const auto spf = build_spf(n + 1);
    u128 ans = 0;

    for (int r = 2; r < n; ++r) {
        const int m = std::min(r - 1, n - r);
        if (m <= 0) {
            continue;
        }

        const auto fa = factorize_small(r - 1, spf);
        const auto fb = factorize_small(r + 1, spf);
        auto factors = merge_factorizations(fa, fb);

        const u64 sq_minus_1 = static_cast<u64>(r) * static_cast<u64>(r) - 1ULL;

        std::function<void(int, u64)> dfs = [&](int idx, u64 cur) {
            if (idx == static_cast<int>(factors.size())) {
                if (cur == 0 || cur > static_cast<u64>(m)) {
                    return;
                }
                const u64 d = cur;
                const u64 p = static_cast<u64>(r) + d;
                const u64 q = static_cast<u64>(r) + sq_minus_1 / d;
                ans += static_cast<u128>(p) + static_cast<u128>(q);
                return;
            }

            const auto [prime, exp] = factors[idx];
            u64 val = cur;
            for (int e = 0; e <= exp; ++e) {
                if (val > static_cast<u64>(m)) {
                    break;
                }
                dfs(idx + 1, val);
                if (e == exp || val > static_cast<u64>(m) / static_cast<u64>(prime)) {
                    break;
                }
                val *= static_cast<u64>(prime);
            }
        };

        dfs(0, 1);
    }

    return ans;
}

static u128 solve_bruteforce(int n) {
    u128 ans = 0;
    for (int p = 2; p <= n; ++p) {
        for (int r = 1; r < p; ++r) {
            const int den = p - r;
            const int num = p * r - 1;
            if (num % den != 0) {
                continue;
            }
            const int q = num / den;
            if (q <= p) {
                continue;
            }
            if ((static_cast<long long>(p) * r) % q != 1) {
                continue;
            }
            if ((static_cast<long long>(q) * r) % p != 1) {
                continue;
            }
            ans += static_cast<u128>(p + q);
        }
    }
    return ans;
}

int main() {
    assert(to_string_u128(solve_fast(5)) == "59");
    assert(to_string_u128(solve_fast(100)) == "697317");
    assert(solve_fast(250) == solve_bruteforce(250));

    const u128 ans = solve_fast(2'000'000);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
