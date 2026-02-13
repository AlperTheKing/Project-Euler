#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

u64 brute_count(int n) {
    std::vector<std::array<std::pair<int, int>, 3>> next(2 * n);
    std::vector<int> deg(2 * n, 0);

    auto id = [&](int size, int ori) {
        return (size - 1) * 2 + ori;
    };

    for (int size = 1; size <= n; ++size) {
        for (int ori = 0; ori < 2; ++ori) {
            const int u = id(size, ori);

            const int nest = (ori == 0 ? size - 1 : size + 1);
            if (1 <= nest && nest <= n) {
                next[u][deg[u]++] = {nest, ori};
            }

            for (int d : {-2, 2}) {
                const int v = size + d;
                if (1 <= v && v <= n) {
                    next[u][deg[u]++] = {v, 1 - ori};
                }
            }
        }
    }

    const u64 full = (u64{1} << n) - 1;
    std::unordered_map<u64, u64> memo;
    memo.reserve(1u << 20);

    std::function<u64(u64, int, int)> dfs = [&](u64 mask, int size, int ori) -> u64 {
        const u64 key = (mask << 6) ^ (u64(size) << 1) ^ u64(ori);
        auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        if (mask == full) {
            return memo[key] = 1;
        }

        const int u = id(size, ori);
        u64 ans = 0;
        for (int i = 0; i < deg[u]; ++i) {
            const int ns = next[u][i].first;
            const int no = next[u][i].second;
            const u64 bit = u64{1} << (ns - 1);
            if ((mask & bit) == 0) {
                ans += dfs(mask | bit, ns, no);
            }
        }

        memo[key] = ans;
        return ans;
    };

    u64 total = 0;
    for (int s = 1; s <= n; ++s) {
        const u64 bit = u64{1} << (s - 1);
        total += dfs(bit, s, 0);
        total += dfs(bit, s, 1);
    }
    return total;
}

i64 S_mod(std::int64_t n) {
    static constexpr i64 MOD = 1'000'000'007LL;

    if (n <= 0) {
        return 0;
    }

    const std::array<i64, 10> base = {0, 2, 2, 6, 12, 16, 22, 36, 58, 82};
    if (n <= 9) {
        return base[static_cast<std::size_t>(n)] % MOD;
    }

    std::array<i64, 10> a = base;
    for (std::int64_t i = 10; i <= n; ++i) {
        i64 v = 0;
        v = (v + 2 * a[9]) % MOD;
        v = (v - 3 * a[8]) % MOD;
        v = (v + 5 * a[7]) % MOD;
        v = (v - 4 * a[6]) % MOD;
        v = (v + 4 * a[5]) % MOD;
        v = (v - 3 * a[4]) % MOD;
        v = (v + a[3]) % MOD;
        v = (v - a[2]) % MOD;
        if (v < 0) {
            v += MOD;
        }

        for (int j = 2; j <= 9; ++j) {
            a[j - 1] = a[j];
        }
        a[9] = v;
    }

    return a[9];
}

void validate() {
    const std::array<u64, 11> exact = {
        0,
        2,
        2,
        6,
        12,
        16,
        22,
        36,
        58,
        82,
        114,
    };

    for (int n = 1; n <= 10; ++n) {
        assert(brute_count(n) == exact[static_cast<std::size_t>(n)]);
        assert(S_mod(n) == static_cast<i64>(exact[static_cast<std::size_t>(n)]));
    }

    assert(S_mod(4) == 12);
    assert(S_mod(8) == 58);
    assert(S_mod(20) == 5560);

    for (int n = 10; n <= 30; ++n) {
        const i64 lhs = S_mod(n);
        i64 rhs = 0;
        rhs += 2 * S_mod(n - 1);
        rhs -= 3 * S_mod(n - 2);
        rhs += 5 * S_mod(n - 3);
        rhs -= 4 * S_mod(n - 4);
        rhs += 4 * S_mod(n - 5);
        rhs -= 3 * S_mod(n - 6);
        rhs += S_mod(n - 7);
        rhs -= S_mod(n - 8);
        rhs %= 1'000'000'007LL;
        if (rhs < 0) {
            rhs += 1'000'000'007LL;
        }
        assert(lhs == rhs);
    }
}

}  // namespace

int main() {
    validate();
    std::cout << S_mod(10'000'000) << '\n';
    return 0;
}
