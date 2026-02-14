#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static std::vector<int> build_spf(int n) {
    std::vector<int> spf(n + 1);
    for (int i = 0; i <= n; ++i) spf[i] = i;
    for (int i = 2; static_cast<long long>(i) * i <= n; ++i) {
        if (spf[i] == i) {
            for (int j = i * i; j <= n; j += i) {
                if (spf[j] == j) spf[j] = i;
            }
        }
    }
    return spf;
}

static void factor_with_spf(int x, const std::vector<int>& spf, std::map<u64, int>& fac) {
    while (x > 1) {
        int p = spf[x];
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        fac[static_cast<u64>(p)] += e;
    }
}

static u128 T_value(int n, const std::vector<int>& spf) {
    int a = n;
    int b = n - 1;

    int t = 0;
    while ((a & 1) == 0) {
        a >>= 1;
        ++t;
    }
    while ((b & 1) == 0) {
        b >>= 1;
        ++t;
    }

    std::map<u64, int> fac;
    factor_with_spf(a, spf, fac);
    factor_with_spf(b, spf, fac);

    std::vector<std::pair<u64, int>> pf(fac.begin(), fac.end());

    u128 ans = 0;
    std::function<void(int, u128)> dfs = [&](int idx, u128 cur) {
        if (idx == static_cast<int>(pf.size())) {
            if (cur > static_cast<u128>(n)) {
                ans += cur - static_cast<u128>(n);
            }

            u128 d2 = cur << t;
            if (d2 > static_cast<u128>(n)) {
                ans += d2 - static_cast<u128>(n);
            }
            return;
        }

        auto [p, e] = pf[idx];
        u128 v = 1;
        for (int k = 0; k <= e; ++k) {
            dfs(idx + 1, cur * v);
            v *= static_cast<u128>(p);
        }
    };

    dfs(0, 1);
    return ans;
}

static u128 U_value(int N) {
    auto spf = build_spf(N);
    u128 total = 0;
    for (int n = 3; n <= N; ++n) {
        total += T_value(n, spf);
    }
    return total;
}

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int d = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + d));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    auto spf100 = build_spf(100);
    assert(T_value(10, spf100) == 148);
    assert(T_value(100, spf100) == 21828);
    assert(U_value(100) == static_cast<u128>(612572));

    const u128 ans = U_value(1'234'567);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
