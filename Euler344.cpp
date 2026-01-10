#include <algorithm>
#include <cstdint>
#include <future>
#include <iostream>
#include <string>
#include <vector>

using u128 = unsigned __int128;

namespace {

constexpr uint64_t kMod = 1000036000099ULL;
constexpr uint64_t kP1 = 1000003ULL;
constexpr uint64_t kP2 = 1000033ULL;

struct SmallComb {
    std::vector<std::vector<unsigned long long>> c;
};

SmallComb build_small_comb(int max_n) {
    SmallComb comb;
    comb.c.assign(max_n + 1, std::vector<unsigned long long>(max_n + 1, 0));
    for (int n = 0; n <= max_n; ++n) {
        comb.c[n][0] = comb.c[n][n] = 1;
        for (int k = 1; k < n; ++k) {
            comb.c[n][k] = comb.c[n - 1][k - 1] + comb.c[n - 1][k];
        }
    }
    return comb;
}

uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    while (exp > 0) {
        if (exp & 1) result = static_cast<uint64_t>((u128)result * base % mod);
        base = static_cast<uint64_t>((u128)base * base % mod);
        exp >>= 1;
    }
    return result;
}

uint64_t mod_inv_prime(uint64_t a, uint64_t p) {
    return mod_pow(a, p - 2, p);
}

struct FactTable {
    std::vector<uint64_t> fact;
    std::vector<uint64_t> invfact;
};

FactTable build_fact_table(int n, uint64_t p) {
    FactTable ft;
    ft.fact.resize(n + 1);
    ft.invfact.resize(n + 1);
    ft.fact[0] = 1;
    for (int i = 1; i <= n; ++i) {
        ft.fact[i] = static_cast<uint64_t>((u128)ft.fact[i - 1] * i % p);
    }
    ft.invfact[n] = mod_pow(ft.fact[n], p - 2, p);
    for (int i = n; i >= 1; --i) {
        ft.invfact[i - 1] = static_cast<uint64_t>((u128)ft.invfact[i] * i % p);
    }
    return ft;
}

uint64_t binom_mod_prime(int n, int k, const FactTable& ft, uint64_t p) {
    if (k < 0 || k > n) return 0;
    uint64_t res = static_cast<uint64_t>((u128)ft.fact[n] * ft.invfact[k] % p);
    res = static_cast<uint64_t>((u128)res * ft.invfact[n - k] % p);
    return res;
}

uint64_t crt(uint64_t a, uint64_t b) {
    uint64_t t = (b + kP2 - (a % kP2)) % kP2;
    uint64_t inv = mod_inv_prime(kP1 % kP2, kP2);
    uint64_t k = static_cast<uint64_t>((u128)t * inv % kP2);
    return static_cast<uint64_t>((u128)a + (u128)kP1 * k);
}

uint64_t binom_mod_composite(int n, int k, const FactTable& f1, const FactTable& f2) {
    uint64_t a = binom_mod_prime(n, k, f1, kP1);
    uint64_t b = binom_mod_prime(n, k, f2, kP2);
    return crt(a, b);
}

struct WaysEvenMod {
    int max_s = 0;
    std::vector<uint64_t> ways;
};

WaysEvenMod build_ways_even_mod(int k, int r, const SmallComb& comb) {
    WaysEvenMod out;
    out.max_s = k + r;
    out.ways.assign(out.max_s + 1, 0);
    for (int x = 0; x <= k; x += 2) {
        unsigned long long cx = comb.c[k][x];
        for (int y = 0; y <= r; ++y) {
            int s = x + y;
            u128 add = static_cast<u128>(cx) * comb.c[r][y];
            uint64_t val = static_cast<uint64_t>(add % kMod);
            uint64_t sum = out.ways[s] + val;
            if (sum >= kMod) sum -= kMod;
            out.ways[s] = sum;
        }
    }
    return out;
}

int bit_length(long long v) {
    int bits = 0;
    while (v > 0) {
        ++bits;
        v >>= 1;
    }
    return bits == 0 ? 1 : bits;
}

uint64_t count_with_ways_mod(long long S, const WaysEvenMod& ways) {
    int max_s = ways.max_s;
    std::vector<uint64_t> dp(max_s + 1, 0);
    std::vector<uint64_t> next(max_s + 1, 0);
    dp[0] = 1;

    int bits = bit_length(S) + 7; // flush remaining carry.
    for (int b = 0; b < bits; ++b) {
        int bit = static_cast<int>((S >> b) & 1LL);
        std::fill(next.begin(), next.end(), 0);
        for (int carry = 0; carry <= max_s; ++carry) {
            uint64_t base = dp[carry];
            if (base == 0) continue;
            for (int s = 0; s <= max_s; ++s) {
                uint64_t ways_s = ways.ways[s];
                if (ways_s == 0) continue;
                int total = s + carry;
                if ((total & 1) != bit) continue;
                int carry_out = total >> 1;
                uint64_t add = static_cast<uint64_t>((u128)base * ways_s % kMod);
                uint64_t sum = next[carry_out] + add;
                if (sum >= kMod) sum -= kMod;
                next[carry_out] = sum;
            }
        }
        dp.swap(next);
    }
    return dp[0];
}

struct WaysEvenExact {
    int max_s = 0;
    std::vector<u128> ways;
};

WaysEvenExact build_ways_even_exact(int k, int r, const SmallComb& comb) {
    WaysEvenExact out;
    out.max_s = k + r;
    out.ways.assign(out.max_s + 1, 0);
    for (int x = 0; x <= k; x += 2) {
        u128 cx = comb.c[k][x];
        for (int y = 0; y <= r; ++y) {
            int s = x + y;
            out.ways[s] += cx * comb.c[r][y];
        }
    }
    return out;
}

u128 count_with_ways_exact(long long S, const WaysEvenExact& ways) {
    int max_s = ways.max_s;
    std::vector<u128> dp(max_s + 1, 0);
    std::vector<u128> next(max_s + 1, 0);
    dp[0] = 1;

    int bits = bit_length(S) + 7;
    for (int b = 0; b < bits; ++b) {
        int bit = static_cast<int>((S >> b) & 1LL);
        std::fill(next.begin(), next.end(), 0);
        for (int carry = 0; carry <= max_s; ++carry) {
            u128 base = dp[carry];
            if (base == 0) continue;
            for (int s = 0; s <= max_s; ++s) {
                u128 ways_s = ways.ways[s];
                if (ways_s == 0) continue;
                int total = s + carry;
                if ((total & 1) != bit) continue;
                int carry_out = total >> 1;
                next[carry_out] += base * ways_s;
            }
        }
        dp.swap(next);
    }
    return dp[0];
}

u128 binom_exact(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    u128 res = 1;
    for (int i = 1; i <= k; ++i) {
        res = res * static_cast<u128>(n - k + i) / static_cast<u128>(i);
    }
    return res;
}

u128 compute_W_exact(int n, int c, const SmallComb& comb) {
    int m = c + 1;
    long long N = n - m;
    u128 total = static_cast<u128>(m) * binom_exact(n, m);
    if (m == 1) return total;

    int k = (m + 1) / 2;
    int r = m - k + 1;
    if (m % 2 == 0) {
        WaysEvenExact ways = build_ways_even_exact(k, r, comb);
        u128 losing = static_cast<u128>(m - 1) * count_with_ways_exact(N, ways);
        return total - losing;
    }

    WaysEvenExact ways_k = build_ways_even_exact(k, r, comb);
    WaysEvenExact ways_km1 = build_ways_even_exact(k - 1, r, comb);
    u128 L0 = count_with_ways_exact(N, ways_k);
    u128 L1 = count_with_ways_exact(N + 1, ways_k) - count_with_ways_exact(N + 1, ways_km1);
    u128 losing = L0 + static_cast<u128>(m - 2) * L1;
    return total - losing;
}

uint64_t compute_W_mod(int n, int c, const SmallComb& comb,
                       const FactTable& f1, const FactTable& f2) {
    int m = c + 1;
    long long N = n - m;
    uint64_t total = static_cast<uint64_t>((u128)m * binom_mod_composite(n, m, f1, f2) % kMod);
    if (m == 1) return total;

    int k = (m + 1) / 2;
    int r = m - k + 1;
    if (m % 2 == 0) {
        WaysEvenMod ways = build_ways_even_mod(k, r, comb);
        uint64_t losing = static_cast<uint64_t>((u128)(m - 1) * count_with_ways_mod(N, ways) % kMod);
        return (total + kMod - losing) % kMod;
    }

    WaysEvenMod ways_k = build_ways_even_mod(k, r, comb);
    WaysEvenMod ways_km1 = build_ways_even_mod(k - 1, r, comb);

    auto fut0 = std::async(std::launch::async, [&]() { return count_with_ways_mod(N, ways_k); });
    auto fut1 = std::async(std::launch::async, [&]() { return count_with_ways_mod(N + 1, ways_k); });
    auto fut2 = std::async(std::launch::async, [&]() { return count_with_ways_mod(N + 1, ways_km1); });

    uint64_t L0 = fut0.get();
    uint64_t L1 = fut1.get();
    uint64_t L1_sub = fut2.get();
    L1 = (L1 + kMod - L1_sub) % kMod;

    uint64_t losing = (L0 + static_cast<uint64_t>((u128)(m - 2) * L1 % kMod)) % kMod;
    return (total + kMod - losing) % kMod;
}

std::string to_string_u128(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v > 0) {
        int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

} // namespace

int main() {
    const int max_small = 52;
    SmallComb comb = build_small_comb(max_small);

    std::cout << "--- Validation Checkpoints ---\n";
    u128 v1 = compute_W_exact(10, 2, comb);
    u128 v2 = compute_W_exact(100, 10, comb);
    std::cout << "W(10, 2) = " << to_string_u128(v1)
              << (v1 == 324 ? " [PASS]" : " [FAIL]") << "\n";
    std::cout << "W(100, 10) = " << to_string_u128(v2)
              << (v2 == 1514704946113500ULL ? " [PASS]" : " [FAIL]") << "\n";

    std::cout << "\n--- Final Solution ---\n";
    const int n = 1000000;
    const int c = 100;
    FactTable f1 = build_fact_table(n, kP1);
    FactTable f2 = build_fact_table(n, kP2);

    uint64_t result = compute_W_mod(n, c, comb, f1, f2);
    std::cout << "W(" << n << ", " << c << ") mod " << kMod << " = " << result << "\n";
    return 0;
}
