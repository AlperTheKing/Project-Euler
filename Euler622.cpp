#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = unsigned long long;
using u128 = __uint128_t;

static u64 mod_mul(u64 a, u64 b, u64 mod) { return (__uint128_t)a * b % mod; }

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    if (mod == 1) return 0;
    u64 r = 1 % mod;
    while (e) {
        if (e & 1) r = mod_mul(r, a, mod);
        a = mod_mul(a, a, mod);
        e >>= 1;
    }
    return r;
}

static bool is_prime(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return n == p;
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](u64 a) -> bool {
        if (a % n == 0) return false;
        u64 x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) return false;
        for (int i = 1; i < s; ++i) {
            x = mod_mul(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };

    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (witness(a)) return false;
    }
    return true;
}

static u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    static std::mt19937_64 rng(987654321);
    std::uniform_int_distribution<u64> dist(0, n - 1);

    while (true) {
        u64 c = dist(rng) % n;
        u64 x = dist(rng) % n;
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) { return (mod_mul(v, v, n) + c) % n; };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            u64 diff = (x > y) ? (x - y) : (y - x);
            d = std::gcd(diff, n);
        }
        if (d != n) return d;
    }
}

static void factor_rec(u64 n, std::map<u64, int> &out) {
    if (n == 1) return;
    if (is_prime(n)) {
        out[n]++;
        return;
    }
    u64 d = pollard_rho(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}

static void gen_divisors(const std::vector<std::pair<u64, int>> &pf, int idx, u64 cur, std::vector<u64> &out) {
    if (idx == (int)pf.size()) {
        out.push_back(cur);
        return;
    }
    auto [p, e] = pf[idx];
    u64 v = 1;
    for (int i = 0; i <= e; ++i) {
        gen_divisors(pf, idx + 1, cur * v, out);
        v *= p;
    }
}

static std::vector<int> prime_factors_distinct(int n) {
    std::vector<int> pf;
    for (int p = 2; 1LL * p * p <= n; ++p) {
        if (n % p != 0) continue;
        pf.push_back(p);
        while (n % p == 0) n /= p;
    }
    if (n > 1) pf.push_back(n);
    return pf;
}

static bool has_exact_order(u64 mod, int ord, const std::vector<int> &ord_primes) {
    if (mod == 1) return false;
    if (mod_pow(2, (u64)ord, mod) != 1) return false;
    for (int p : ord_primes) {
        if (mod_pow(2, (u64)(ord / p), mod) == 1) return false;
    }
    return true;
}

static u128 sum_n_with_shuffle_order(int ord) {
    const u64 M = (ord == 64) ? ~0ULL : ((1ULL << ord) - 1ULL);
    std::map<u64, int> f;
    factor_rec(M, f);

    std::vector<std::pair<u64, int>> pf;
    pf.reserve(f.size());
    for (auto [p, e] : f) pf.push_back({p, e});

    std::vector<u64> divs;
    divs.reserve(8192);
    gen_divisors(pf, 0, 1, divs);

    const std::vector<int> ord_primes = prime_factors_distinct(ord);
    u128 sum = 0;
    for (u64 d : divs) {
        if (has_exact_order(d, ord, ord_primes)) sum += (u128)(d + 1);
    }
    return sum;
}

static std::string to_string_u128(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) {
        int digit = (int)(x % 10);
        s.push_back(char('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

int main() {
    assert(sum_n_with_shuffle_order(8) == 412);
    std::cout << to_string_u128(sum_n_with_shuffle_order(60)) << "\n";
    return 0;
}

