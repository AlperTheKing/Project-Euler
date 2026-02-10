#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

static constexpr u64 kLimit = 10'000'000'000'000'000ULL;  // 1e16

static u64 gcd_u64(u64 a, u64 b) {
    while (b) {
        u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static u64 mul_mod_u64(u64 a, u64 b, u64 mod) { return static_cast<u64>((static_cast<u128>(a) * b) % mod); }

static u64 pow_mod_u64(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) {
        if (d & 1) r = mul_mod_u64(r, a, mod);
        a = mul_mod_u64(a, a, mod);
        d >>= 1;
    }
    return r;
}

static bool is_prime_u64(u64 n) {
    if (n < 2) return false;
    static u64 small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }

    u64 d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](u64 a) -> bool {
        if (a % n == 0) return false;
        u64 x = pow_mod_u64(a, d, n);
        if (x == 1 || x == n - 1) return false;
        for (u64 i = 1; i < s; ++i) {
            x = mul_mod_u64(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };

    static u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (witness(a)) return false;
    }
    return true;
}

static u64 splitmix64(u64& x) {
    u64 z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static thread_local u64 rng_state = 0x123456789abcdef0ULL;

static u64 rand_u64(u64 lo, u64 hi) {
    u64 r = splitmix64(rng_state);
    return lo + (hi > lo ? (r % (hi - lo + 1)) : 0);
}

static u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    if (n % 3ULL == 0) return 3;

    while (true) {
        u64 c = rand_u64(1, n - 1);
        u64 x = rand_u64(0, n - 1);
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) { return (mul_mod_u64(v, v, n) + c) % n; };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            u64 diff = (x > y) ? (x - y) : (y - x);
            d = gcd_u64(diff, n);
        }
        if (d != n) return d;
    }
}

static void factor_rec(u64 n, std::vector<u64>& out) {
    if (n == 1) return;
    if (is_prime_u64(n)) {
        out.push_back(n);
        return;
    }
    u64 d = pollard_rho(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}

static std::vector<int> exponent_pattern(u64 n) {
    if (n == 1) return {};
    std::vector<u64> fac;
    factor_rec(n, fac);
    std::sort(fac.begin(), fac.end());
    std::vector<int> exps;
    for (std::size_t i = 0; i < fac.size();) {
        std::size_t j = i;
        while (j < fac.size() && fac[j] == fac[i]) ++j;
        exps.push_back(static_cast<int>(j - i));
        i = j;
    }
    std::sort(exps.begin(), exps.end(), std::greater<int>());
    return exps;
}

static u128 binom_u128(int n, int k) {
    if (k < 0 || k > n) return 0;
    k = std::min(k, n - k);
    u128 r = 1;
    for (int i = 1; i <= k; ++i) {
        r = (r * static_cast<u128>(n - k + i)) / static_cast<u128>(i);
    }
    return r;
}

static u64 g_from_pattern(const std::vector<int>& exps, u64 limit) {
    if (exps.empty()) return 1;
    int omega = 0;
    for (int e : exps) omega += e;

    std::vector<std::vector<u128>> C(omega + 1, std::vector<u128>(omega + 1, 0));
    C[0][0] = 1;
    for (int n = 1; n <= omega; ++n) {
        C[n][0] = C[n][n] = 1;
        for (int k = 1; k < n; ++k) C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
    }

    std::vector<u128> total(static_cast<std::size_t>(omega + 1), 0);
    for (int k = 1; k <= omega; ++k) {
        u128 prod = 1;
        for (int a : exps) {
            prod *= binom_u128(a + k - 1, k - 1);
        }
        total[static_cast<std::size_t>(k)] = prod;
    }

    u128 g = 0;
    for (int k = 1; k <= omega; ++k) {
        i128 fk = 0;
        for (int i = 1; i <= k; ++i) {
            const u128 term = C[k][i] * total[static_cast<std::size_t>(i)];
            if (((k - i) & 1) != 0) fk -= static_cast<i128>(term);
            else fk += static_cast<i128>(term);
        }
        g += static_cast<u128>(fk);
        if (g > static_cast<u128>(limit)) return limit + 1;
    }
    return static_cast<u64>(g);
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

struct Enumerator {
    std::vector<u64> primes;
    std::vector<int> cur;
    std::unordered_map<u64, std::vector<int>> pattern_cache;
    std::vector<u64> solutions;

    Enumerator() {
        // Enough primes for any exponent pattern with minimal value <= 1e16.
        primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
        pattern_cache.reserve(1 << 16);
        solutions.reserve(64);
    }

    const std::vector<int>& cached_pattern(u64 n) {
        auto it = pattern_cache.find(n);
        if (it != pattern_cache.end()) return it->second;
        std::vector<int> p = exponent_pattern(n);
        auto [ins_it, _] = pattern_cache.emplace(n, std::move(p));
        return ins_it->second;
    }

    void dfs(int last_exp, u128 min_n) {
        const u64 g = g_from_pattern(cur, kLimit);
        if (g > kLimit) return;

        if (static_cast<u128>(g) >= min_n) {
            const std::vector<int>& pat = cached_pattern(g);
            if (pat == cur) solutions.push_back(g);
        }

        const std::size_t idx = cur.size();
        if (idx >= primes.size()) return;
        const u64 p = primes[idx];

        u128 pow = 1;
        for (int e = 1; e <= last_exp; ++e) {
            pow *= static_cast<u128>(p);
            const u128 next_min = min_n * pow;
            if (next_min > kLimit) break;
            cur.push_back(e);
            dfs(e, next_min);
            cur.pop_back();
        }
    }

    std::vector<u64> run() {
        cur.clear();
        dfs(53, 1);
        std::sort(solutions.begin(), solutions.end());
        solutions.erase(std::unique(solutions.begin(), solutions.end()), solutions.end());
        return solutions;
    }
};

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Statement checkpoints.
    assert(g_from_pattern(exponent_pattern(12), kLimit) == 8ULL);
    assert(g_from_pattern(exponent_pattern(48), kLimit) == 48ULL);
    assert(g_from_pattern(exponent_pattern(120), kLimit) == 132ULL);

    Enumerator en;
    const std::vector<u64> sols = en.run();

    u128 sum = 0;
    for (u64 x : sols) sum += static_cast<u128>(x);
    std::cout << to_string_u128(sum) << '\n';
    return 0;
}

