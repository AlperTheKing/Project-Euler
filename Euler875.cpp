#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr int LIMIT = 12'345'678;
static constexpr u64 MOD = 1'001'961'001ULL;

static u64 q_prime_power_mod(int p, int e) {
    if (e == 0) return 1;

    if (p == 2) {
        u64 q = 128 % MOD;
        u64 add = 2'048 % MOD;
        for (int k = 2; k <= e; ++k) {
            q = (q * 128ULL + add) % MOD;
            add = (add * 16ULL) % MOD;
        }
        return q;
    }

    u64 pm = static_cast<u64>(p) % MOD;
    u64 p2 = (pm * pm) % MOD;
    u64 p3 = (p2 * pm) % MOD;
    u64 p4 = (p2 * p2) % MOD;
    u64 p7 = (p4 * p3) % MOD;
    u64 coeff = (pm + MOD - 1) % MOD;

    u64 q = 1;
    u64 term = p3;
    for (int k = 1; k <= e; ++k) {
        q = (q * p7 + coeff * term) % MOD;
        term = (term * p4) % MOD;
    }
    return q;
}

static u128 q_prime_power_exact(int p, int e) {
    if (e == 0) return 1;

    if (p == 2) {
        u128 q = 128;
        u128 add = 2048;
        for (int k = 2; k <= e; ++k) {
            q = q * 128 + add;
            add *= 16;
        }
        return q;
    }

    u128 p7 = 1;
    for (int i = 0; i < 7; ++i) p7 *= static_cast<u128>(p);
    u128 p4 = 1;
    for (int i = 0; i < 4; ++i) p4 *= static_cast<u128>(p);
    u128 term = static_cast<u128>(p) * static_cast<u128>(p) * static_cast<u128>(p);

    u128 q = 1;
    for (int k = 1; k <= e; ++k) {
        q = q * p7 + static_cast<u128>(p - 1) * term;
        term *= p4;
    }
    return q;
}

static u128 q_exact_from_factorization(int n) {
    int x = n;
    u128 q = 1;
    for (int p = 2; 1LL * p * p <= x; ++p) {
        if (x % p != 0) continue;
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        q *= q_prime_power_exact(p, e);
    }
    if (x > 1) q *= q_prime_power_exact(x, 1);
    return q;
}

static u64 q_bruteforce_small(int n) {
    std::vector<int> sq(n);
    for (int i = 0; i < n; ++i) sq[i] = static_cast<int>((1LL * i * i) % n);

    std::vector<u64> cnt(n, 0);
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
            for (int c = 0; c < n; ++c) {
                for (int d = 0; d < n; ++d) {
                    int s = sq[a] + sq[b] + sq[c] + sq[d];
                    s %= n;
                    ++cnt[s];
                }
            }
        }
    }

    u64 q = 0;
    for (u64 v : cnt) q += v * v;
    return q;
}

static u64 compute_Q_mod(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            long long v = 1LL * p * i;
            if (v > n) break;
            spf[static_cast<int>(v)] = p;
            if (p == spf[i]) break;
        }
    }

    u64 sum = 1;
    for (int i = 2; i <= n; ++i) {
        int x = i;
        u64 qn = 1;
        while (x > 1) {
            int p = spf[x];
            int e = 0;
            do {
                x /= p;
                ++e;
            } while (x > 1 && spf[x] == p);

            qn = (qn * q_prime_power_mod(p, e)) % MOD;
        }

        sum += qn;
        if (sum >= MOD) sum -= MOD;
    }

    return sum;
}

int main() {
    assert(q_bruteforce_small(4) == 18'432ULL);
    for (int n = 1; n <= 8; ++n) {
        assert(q_bruteforce_small(n) == static_cast<u64>(q_exact_from_factorization(n)));
    }

    u128 q10 = 0;
    for (int i = 1; i <= 10; ++i) q10 += q_exact_from_factorization(i);
    assert(static_cast<u64>(q10) == 18'573'381ULL);
    assert(compute_Q_mod(10) == 18'573'381ULL);

    std::cout << compute_Q_mod(LIMIT) << '\n';
    return 0;
}
