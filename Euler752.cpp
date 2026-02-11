#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Pair {
    u64 a;
    u64 b;
};

Pair mul(const Pair x, const Pair y, const u64 mod) {
    return {
        (x.a * y.a + 7ULL * x.b * y.b) % mod,
        (x.a * y.b + x.b * y.a) % mod,
    };
}

Pair pow_u(const u64 exp, const u64 mod) {
    Pair result{1ULL, 0ULL};
    Pair base{1ULL, 1ULL};
    u64 e = exp;
    while (e > 0ULL) {
        if ((e & 1ULL) != 0ULL) {
            result = mul(result, base, mod);
        }
        base = mul(base, base, mod);
        e >>= 1ULL;
    }
    return result;
}

std::vector<int> build_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1));
    for (int i = 0; i <= n; ++i) {
        spf[static_cast<std::size_t>(i)] = i;
    }
    for (int i = 2; i * i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            if (spf[static_cast<std::size_t>(j)] == j) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return spf;
}

std::vector<int> distinct_prime_factors(int n, const std::vector<int>& spf) {
    std::vector<int> factors;
    while (n > 1) {
        const int p = spf[static_cast<std::size_t>(n)];
        factors.push_back(p);
        while (n % p == 0) {
            n /= p;
        }
    }
    return factors;
}

u64 order_mod_prime(const int p, const std::vector<int>& spf) {
    if (p == 7) {
        return 7ULL;
    }

    const u64 ls = [&]() -> u64 {
        u64 r = 1ULL;
        u64 b = 7ULL % static_cast<u64>(p);
        u64 e = static_cast<u64>((p - 1) / 2);
        const u64 mod = static_cast<u64>(p);
        while (e > 0ULL) {
            if ((e & 1ULL) != 0ULL) {
                r = (r * b) % mod;
            }
            b = (b * b) % mod;
            e >>= 1ULL;
        }
        return r;
    }();

    u64 exponent = 0ULL;
    std::vector<int> factors;
    if (ls == 1ULL) {
        exponent = static_cast<u64>(p - 1);
        factors = distinct_prime_factors(p - 1, spf);
    } else {
        exponent = static_cast<u64>(p - 1) * static_cast<u64>(p + 1);
        factors = distinct_prime_factors(p - 1, spf);
        std::vector<int> extra = distinct_prime_factors(p + 1, spf);
        factors.insert(factors.end(), extra.begin(), extra.end());
        std::sort(factors.begin(), factors.end());
        factors.erase(std::unique(factors.begin(), factors.end()), factors.end());
    }

    u64 ord = exponent;
    for (const int q : factors) {
        const u64 qq = static_cast<u64>(q);
        while (ord % qq == 0ULL) {
            const Pair test = pow_u(ord / qq, static_cast<u64>(p));
            if (test.a == 1ULL && test.b == 0ULL) {
                ord /= qq;
            } else {
                break;
            }
        }
    }
    return ord;
}

u64 G(const int n) {
    const int limit = n + 1;
    const std::vector<int> spf = build_spf(limit);

    std::vector<u64> ord_prime_power(static_cast<std::size_t>(n + 1), 0ULL);

    for (int p = 5; p <= n; ++p) {
        if (spf[static_cast<std::size_t>(p)] != p) {
            continue;
        }

        u64 ord = order_mod_prime(p, spf);
        int pk = p;
        ord_prime_power[static_cast<std::size_t>(pk)] = ord;

        while (pk <= n / p) {
            pk *= p;
            const Pair test = pow_u(ord, static_cast<u64>(pk));
            if (!(test.a == 1ULL && test.b == 0ULL)) {
                ord *= static_cast<u64>(p);
            }
            ord_prime_power[static_cast<std::size_t>(pk)] = ord;
        }
    }

    u64 total = 0ULL;
    for (int x = 2; x <= n; ++x) {
        if ((x % 2) == 0 || (x % 3) == 0) {
            continue;
        }

        int t = x;
        u64 gx = 1ULL;
        while (t > 1) {
            const int p = spf[static_cast<std::size_t>(t)];
            int pk = 1;
            while (t % p == 0) {
                t /= p;
                pk *= p;
            }
            gx = std::lcm(gx, ord_prime_power[static_cast<std::size_t>(pk)]);
        }
        total += gx;
    }

    return total;
}

}  // namespace

int main() {
    assert(G(100) == 28'891ULL);
    assert(G(1'000) == 13'131'583ULL);

    std::cout << G(1'000'000) << '\n';
    return 0;
}
