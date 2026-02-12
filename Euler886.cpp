#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

static constexpr int MOD = 83'456'729;

static int odd_part(int x) {
    while ((x & 1) == 0) x >>= 1;
    return x;
}

static int P_count(int n) {
    std::vector<int> evens;
    std::vector<int> odds;
    for (int x = 2; x <= n; ++x) {
        if (x & 1) odds.push_back(x);
        else evens.push_back(x);
    }

    std::vector<std::pair<int, int>> even_types;
    for (int e : evens) {
        int d = odd_part(e);
        auto it = std::find_if(even_types.begin(), even_types.end(), [&](const auto& p) { return p.first == d; });
        if (it == even_types.end()) even_types.push_back({d, 1});
        else ++it->second;
    }
    std::sort(even_types.begin(), even_types.end());

    int E = static_cast<int>(even_types.size());
    std::vector<int> e_part(E), e_count(E);
    for (int i = 0; i < E; ++i) {
        e_part[i] = even_types[i].first;
        e_count[i] = even_types[i].second;
    }

    std::vector<u64> o_masks;
    std::vector<int> o_count;
    for (int o : odds) {
        u64 mask = 0;
        for (int i = 0; i < E; ++i) {
            if (std::gcd(o, e_part[i]) == 1) mask |= (1ULL << i);
        }

        auto it = std::find(o_masks.begin(), o_masks.end(), mask);
        if (it == o_masks.end()) {
            o_masks.push_back(mask);
            o_count.push_back(1);
        } else {
            int idx = static_cast<int>(it - o_masks.begin());
            ++o_count[idx];
        }
    }

    int O = static_cast<int>(o_masks.size());

    std::vector<int> baseE(E), baseO(O);
    for (int i = 0; i < E; ++i) baseE[i] = e_count[i] + 1;
    for (int j = 0; j < O; ++j) baseO[j] = o_count[j] + 1;

    std::vector<u64> radixE(E, 1), radixO(O, 1);
    for (int i = 1; i < E; ++i) radixE[i] = radixE[i - 1] * static_cast<u64>(baseE[i - 1]);
    for (int j = 1; j < O; ++j) radixO[j] = radixO[j - 1] * static_cast<u64>(baseO[j - 1]);

    u64 totalOStates = 1;
    for (int j = 0; j < O; ++j) totalOStates *= static_cast<u64>(baseO[j]);

    auto getE = [&](u64 enc, int i) -> int {
        return static_cast<int>((enc / radixE[i]) % static_cast<u64>(baseE[i]));
    };
    auto getO = [&](u64 enc, int j) -> int {
        return static_cast<int>((enc / radixO[j]) % static_cast<u64>(baseO[j]));
    };

    u64 encE0 = 0, encO0 = 0;
    for (int i = 0; i < E; ++i) encE0 += static_cast<u64>(e_count[i]) * radixE[i];
    for (int j = 0; j < O; ++j) encO0 += static_cast<u64>(o_count[j]) * radixO[j];

    std::unordered_map<u64, int> memo;
    memo.reserve(2'000'000);

    std::function<int(u64, u64, int)> dp = [&](u64 encE, u64 encO, int last) -> int {
        if (encO == 0) return 1;

        u64 key = ((encE * totalOStates + encO) * static_cast<u64>(E)) + static_cast<u64>(last);
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;

        long long ans = 0;

        for (int j = 0; j < O; ++j) {
            int cj = getO(encO, j);
            if (cj == 0) continue;
            if (((o_masks[j] >> last) & 1ULL) == 0) continue;

            u64 encO2 = encO - radixO[j];

            for (int i = 0; i < E; ++i) {
                int ci = getE(encE, i);
                if (ci == 0) continue;
                if (((o_masks[j] >> i) & 1ULL) == 0) continue;

                u64 encE2 = encE - radixE[i];
                long long ways = (static_cast<long long>(cj) * ci) % MOD;
                ans += ways * dp(encE2, encO2, i);
                ans %= MOD;
            }
        }

        int ret = static_cast<int>(ans % MOD);
        memo[key] = ret;
        return ret;
    };

    long long total = 0;
    for (int i = 0; i < E; ++i) {
        int ci = getE(encE0, i);
        if (ci == 0) continue;
        u64 encE1 = encE0 - radixE[i];
        total += (static_cast<long long>(ci) * dp(encE1, encO0, i)) % MOD;
        total %= MOD;
    }

    return static_cast<int>(total % MOD);
}

static int brute_count(int n) {
    std::vector<int> v;
    for (int x = 2; x <= n; ++x) v.push_back(x);

    int cnt = 0;
    std::sort(v.begin(), v.end());
    do {
        bool ok = true;
        for (int i = 1; i < static_cast<int>(v.size()); ++i) {
            if (std::gcd(v[i - 1], v[i]) != 1) {
                ok = false;
                break;
            }
        }
        if (ok) ++cnt;
    } while (std::next_permutation(v.begin(), v.end()));

    return cnt;
}

int main() {
    assert(P_count(4) == 2);
    assert(P_count(10) == 576);
    assert(P_count(6) == brute_count(6));
    assert(P_count(8) == brute_count(8));

    std::cout << P_count(34) << '\n';
    return 0;
}
