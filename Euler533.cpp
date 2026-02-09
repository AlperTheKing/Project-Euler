#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kLast = 1'000'000'000ULL;  // last 9 digits

u64 pow_mod(u64 a, int e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1) {
            r = static_cast<u128>(r) * a % mod;
        }
        a = static_cast<u128>(a) * a % mod;
        e >>= 1;
    }
    return r;
}

std::vector<bool> sieve_is_prime(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; (u64)p * p <= static_cast<u64>(n); ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int m = p * p; m <= n; m += p) {
            is_prime[static_cast<std::size_t>(m)] = false;
        }
    }
    return is_prime;
}

struct BestState {
    long double log_value = -1e300L;
    u64 mod_value = 0;
    u64 L = 0;
};

u64 exponent_in_factorization(const std::vector<std::pair<u64, int>>& fac, const u64 p) {
    for (const auto& [q, e] : fac) {
        if (q == p) {
            return static_cast<u64>(e);
        }
    }
    return 0;
}

void evaluate_L(const u64 M, const std::vector<bool>& is_prime, const u64 L,
                const std::vector<std::pair<u64, int>>& fac, BestState& best) {
    // Build full divisor list of L.
    std::vector<u64> divisors;
    divisors.reserve(4096);
    divisors.push_back(1);
    for (const auto& [p, e] : fac) {
        const std::size_t cur = divisors.size();
        u64 pe = 1;
        for (int i = 1; i <= e; ++i) {
            pe *= p;
            for (std::size_t j = 0; j < cur; ++j) {
                divisors.push_back(divisors[j] * pe);
            }
        }
    }

    const u64 v2 = exponent_in_factorization(fac, 2);
    const int exp2 = (v2 == 0) ? 1 : static_cast<int>(v2 + 2);

    long double logN = static_cast<long double>(exp2) * std::log(2.0L);
    u64 modN = pow_mod(2, exp2, kLast);

    for (const u64 d : divisors) {
        const u64 p = d + 1;
        if (p == 2) {
            continue;
        }
        if (p > M + 1) {
            continue;
        }
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        const int exp = static_cast<int>(exponent_in_factorization(fac, p)) + 1;
        logN += static_cast<long double>(exp) * std::log(static_cast<long double>(p));
        modN = static_cast<u128>(modN) * pow_mod(p, exp, kLast) % kLast;
    }

    if (logN > best.log_value) {
        best.log_value = logN;
        best.mod_value = modN;
        best.L = L;
    }
}

void dfs_generate(const std::vector<u64>& primes, const u64 M, const std::vector<bool>& is_prime, int idx, int max_exp,
                  u64 cur, std::vector<std::pair<u64, int>>& fac, BestState& best) {
    evaluate_L(M, is_prime, cur, fac, best);
    if (idx >= static_cast<int>(primes.size())) {
        return;
    }

    const u64 p = primes[static_cast<std::size_t>(idx)];
    u64 val = cur;
    for (int e = 1; e <= max_exp; ++e) {
        if (val > M / p) {
            break;
        }
        val *= p;
        fac.push_back({p, e});
        dfs_generate(primes, M, is_prime, idx + 1, e, val, fac, best);
        fac.pop_back();
    }
}

u64 solve_last9(const u64 bound_minus_1) {
    const u64 M = bound_minus_1;
    const auto is_prime = sieve_is_prime(static_cast<int>(M + 1));
    // Enough primes to generate all relevant "highly composite"-style candidates under 2e7.
    const std::vector<u64> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43};

    BestState best;
    std::vector<std::pair<u64, int>> fac;
    dfs_generate(primes, M, is_prime, 0, 60, 1, fac, best);

    // L(n) = max{ m : lambda(m) < n } + 1, and last 9 digits are requested.
    return (best.mod_value + 1) % kLast;
}

bool run_checkpoints() {
    if (solve_last9(6 - 1) != 241ULL) {
        std::cerr << "Checkpoint failed: L(6)\n";
        return false;
    }
    if (solve_last9(100 - 1) != 174'525'281ULL) {
        std::cerr << "Checkpoint failed: L(100) last 9 digits\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }
    const u64 last9 = solve_last9(20'000'000ULL - 1);
    std::cout << std::setw(9) << std::setfill('0') << last9 << '\n';
    return 0;
}

