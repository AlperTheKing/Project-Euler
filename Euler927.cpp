#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

u64 mod_pow(u64 base, int exp, int mod) {
    if (mod == 1) {
        return 0;
    }
    u64 result = 1 % mod;
    base %= static_cast<u64>(mod);
    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % static_cast<u64>(mod);
        }
        base = (base * base) % static_cast<u64>(mod);
        exp >>= 1;
    }
    return result;
}

u64 next_term(u64 x, int p, int mod) {
    if (p == 2) {
        return (x * x + 1ULL) % static_cast<u64>(mod);
    }
    return (mod_pow(x, p, mod) + 1ULL) % static_cast<u64>(mod);
}

bool hits_zero(int mod, int p) {
    if (mod == 1) {
        return true;
    }

    u64 tortoise = 1ULL % static_cast<u64>(mod);
    u64 hare = tortoise;
    if (tortoise == 0) {
        return true;
    }

    while (true) {
        tortoise = next_term(tortoise, p, mod);
        if (tortoise == 0) {
            return true;
        }

        hare = next_term(hare, p, mod);
        if (hare == 0) {
            return true;
        }

        hare = next_term(hare, p, mod);
        if (hare == 0) {
            return true;
        }

        if (tortoise == hare) {
            break;
        }
    }

    u64 x = next_term(tortoise, p, mod);
    while (x != tortoise) {
        if (x == 0) {
            return true;
        }
        x = next_term(x, p, mod);
    }

    return false;
}

struct Sieve {
    std::vector<int> primes;
    std::vector<int> spf;
};

Sieve build_sieve(int n) {
    std::vector<int> spf(n + 1, 0);
    std::vector<int> primes;
    primes.reserve(n / 10);

    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(p) * i;
            if (v > n || p > spf[i]) {
                break;
            }
            spf[static_cast<int>(v)] = p;
        }
    }

    return {std::move(primes), std::move(spf)};
}

std::vector<int> unique_prime_divisors(int x, const std::vector<int>& spf) {
    std::vector<int> divisors;
    while (x > 1) {
        const int p = spf[x];
        divisors.push_back(p);
        while (x % p == 0) {
            x /= p;
        }
    }
    return divisors;
}

std::vector<int> good_primes_upto(int limit) {
    const Sieve sieve = build_sieve(limit);

    std::vector<int> good;
    for (int q : sieve.primes) {
        if (q == 2) {
            good.push_back(q);
            continue;
        }

        if (!hits_zero(q, 2)) {
            continue;
        }

        bool ok = true;
        const std::vector<int> divisors = unique_prime_divisors(q - 1, sieve.spf);
        for (int p : divisors) {
            if (p == 2) {
                continue;
            }
            if (!hits_zero(q, p)) {
                ok = false;
                break;
            }
        }

        if (ok) {
            good.push_back(q);
        }
    }

    return good;
}

u64 sum_squarefree_products_leq(const std::vector<int>& primes, i64 limit) {
    u64 total = 0;

    auto dfs = [&](auto&& self, std::size_t idx, i64 current) -> void {
        total += static_cast<u64>(current);
        for (std::size_t i = idx; i < primes.size(); ++i) {
            const i64 p = primes[i];
            if (current > limit / p) {
                break;
            }
            self(self, i + 1, current * p);
        }
    };

    dfs(dfs, 0, 1);
    return total;
}

u64 solve(i64 n) {
    const std::vector<int> good = good_primes_upto(static_cast<int>(n));

    for (int q : good) {
        if (static_cast<i64>(q) * q > n) {
            continue;
        }
        assert(!hits_zero(q * q, 2));
    }

    return sum_squarefree_products_leq(good, n);
}

void run_validations() {
    assert(solve(20) == 18);
    assert(solve(1000) == 2089);
}

}  // namespace

int main() {
    run_validations();
    constexpr i64 kN = 10'000'000;
    std::cout << solve(kN) << '\n';
    return 0;
}
