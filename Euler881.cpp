#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static int width_from_exponents(const std::vector<int>& exps) {
    std::vector<int> poly(1, 1);
    for (int e : exps) {
        std::vector<int> next(static_cast<std::size_t>(poly.size() + e), 0);
        for (std::size_t i = 0; i < poly.size(); ++i) {
            int v = poly[i];
            for (int t = 0; t <= e; ++t) {
                next[i + static_cast<std::size_t>(t)] += v;
            }
        }
        poly.swap(next);
    }
    return *std::max_element(poly.begin(), poly.end());
}

static std::vector<int> factor_exponents(u64 n) {
    std::vector<int> exps;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) continue;
        int e = 0;
        while (n % p == 0) {
            n /= p;
            ++e;
        }
        exps.push_back(e);
    }
    if (n > 1) exps.push_back(1);
    std::sort(exps.begin(), exps.end(), std::greater<int>());
    return exps;
}

static int g_of_n(u64 n) {
    return width_from_exponents(factor_exponents(n));
}

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct Solver {
    static constexpr int TARGET = 10'000;
    std::vector<int> primes;
    u128 best;

    Solver() : primes{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61}, best(0) {
        best = 1;
        for (int i = 0; i < 16; ++i) best *= static_cast<u128>(primes[i]);
    }

    void dfs(int pos, int max_exp, u128 n, std::vector<int>& exps) {
        int w = width_from_exponents(exps);
        if (w >= TARGET) {
            if (n < best) best = n;
            return;
        }
        if (pos >= static_cast<int>(primes.size())) return;

        u128 mult = 1;
        u128 p = static_cast<u128>(primes[pos]);
        for (int e = 1; e <= max_exp; ++e) {
            mult *= p;
            u128 n2 = n * mult;
            if (n2 >= best) break;
            exps.push_back(e);
            dfs(pos + 1, e, n2, exps);
            exps.pop_back();
        }
    }

    u128 solve() {
        std::vector<int> exps;
        dfs(0, 64, 1, exps);
        return best;
    }
};

int main() {
    assert(g_of_n(45) == 2);
    assert(g_of_n(5040) == 12);

    Solver solver;
    u128 ans = solver.solve();
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
