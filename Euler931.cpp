#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr u64 kMod = 715'827'883ULL;

i128 tri(u64 x) {
    return static_cast<i128>(x) * static_cast<i128>(x + 1) / 2;
}

bool is_prime_small(int x) {
    if (x < 2) {
        return false;
    }
    for (int d = 2; static_cast<i64>(d) * d <= x; ++d) {
        if (x % d == 0) {
            return false;
        }
    }
    return true;
}

int phi_small(int n) {
    int result = n;
    int x = n;
    for (int p = 2; static_cast<i64>(p) * p <= x; ++p) {
        if (x % p != 0) {
            continue;
        }
        while (x % p == 0) {
            x /= p;
        }
        result -= result / p;
    }
    if (x > 1) {
        result -= result / x;
    }
    return result;
}

i64 t_bruteforce_graph(int n) {
    std::vector<int> divisors;
    for (int d = 1; static_cast<i64>(d) * d <= n; ++d) {
        if (n % d != 0) {
            continue;
        }
        divisors.push_back(d);
        if (d * d != n) {
            divisors.push_back(n / d);
        }
    }

    i64 total = 0;
    for (int a : divisors) {
        for (int b : divisors) {
            if (a % b != 0) {
                continue;
            }
            const int q = a / b;
            if (!is_prime_small(q)) {
                continue;
            }
            total += static_cast<i64>(phi_small(a) - phi_small(b));
        }
    }
    return total;
}

i64 t_formula_small(int n) {
    int x = n;
    i64 total = 0;

    for (int p = 2; static_cast<i64>(p) * p <= x; ++p) {
        if (x % p != 0) {
            continue;
        }
        i64 pe = 1;
        while (x % p == 0) {
            x /= p;
            pe *= p;
        }
        total += static_cast<i64>(n) - n / p - n / pe;
    }
    if (x > 1) {
        total += static_cast<i64>(n) - 2LL * (n / x);
    }

    return total;
}

struct PrimePrefix {
    u64 n = 0;
    u64 root = 0;

    std::vector<u64> values;
    std::vector<int> idx_small;
    std::vector<int> idx_large;

    std::vector<i64> g_count;
    std::vector<i128> g_sum;

    std::vector<int> primes;

    int id(u64 x) const {
        if (x <= root) {
            return idx_small[x];
        }
        return idx_large[n / x];
    }

    explicit PrimePrefix(u64 limit) : n(limit), root(static_cast<u64>(std::sqrt(static_cast<long double>(limit)))) {
        for (u64 l = 1, r; l <= n; l = r + 1) {
            const u64 v = n / l;
            r = n / v;
            values.push_back(v);
        }

        idx_small.assign(root + 1, -1);
        idx_large.assign(root + 1, -1);

        const int m = static_cast<int>(values.size());
        g_count.resize(m);
        g_sum.resize(m);

        for (int i = 0; i < m; ++i) {
            const u64 v = values[i];
            if (v <= root) {
                idx_small[v] = i;
            } else {
                idx_large[n / v] = i;
            }
            g_count[i] = static_cast<i64>(v) - 1;
            g_sum[i] = tri(v) - 1;
        }

        std::vector<int> lp(root + 1, 0);
        for (u64 i = 2; i <= root; ++i) {
            if (lp[i] == 0) {
                lp[i] = static_cast<int>(i);
                primes.push_back(static_cast<int>(i));
            }
            for (int p : primes) {
                const u64 v = i * static_cast<u64>(p);
                if (v > root || p > lp[i]) {
                    break;
                }
                lp[v] = p;
            }
        }

        for (int p : primes) {
            const u64 p64 = static_cast<u64>(p);
            const u64 p2 = p64 * p64;
            if (p2 > n) {
                break;
            }

            const int idx_pm1 = id(p64 - 1);

            int i = 0;
            while (i < m && values[i] >= p2) {
                const int j = id(values[i] / p64);
                g_count[i] -= g_count[j] - g_count[idx_pm1];
                g_sum[i] -= static_cast<i128>(p64) * (g_sum[j] - g_sum[idx_pm1]);
                ++i;
            }
        }
    }

    i64 prime_count(u64 x) const {
        if (x < 2) {
            return 0;
        }
        return g_count[id(x)];
    }

    i128 prime_sum(u64 x) const {
        if (x < 2) {
            return 0;
        }
        return g_sum[id(x)];
    }
};

u64 solve(u64 N) {
    PrimePrefix prefix(N);
    const u64 root = prefix.root;

    u64 answer = 0;

    auto add_mod_i128 = [&](i128 x) {
        x %= static_cast<i128>(kMod);
        if (x < 0) {
            x += static_cast<i128>(kMod);
        }
        answer += static_cast<u64>(x);
        if (answer >= kMod) {
            answer -= kMod;
        }
    };

    for (int p : prefix.primes) {
        const u64 p64 = static_cast<u64>(p);
        u64 pe = p64;

        while (pe <= N) {
            const u64 x = N / pe;
            const i128 s = tri(x) - static_cast<i128>(p64) * tri(x / p64);
            const i128 coef = static_cast<i128>(pe) - static_cast<i128>(pe / p64) - 1;
            add_mod_i128(coef * s);

            if (pe > N / p64) {
                break;
            }
            pe *= p64;
        }
    }

    const u64 max_q = N / (root + 1);
    for (u64 q = 1; q <= max_q; ++q) {
        u64 L = N / (q + 1) + 1;
        u64 R = N / q;

        if (R <= root) {
            break;
        }
        if (L <= root) {
            L = root + 1;
        }
        if (L > R) {
            continue;
        }

        const i64 cnt = prefix.prime_count(R) - prefix.prime_count(L - 1);
        if (cnt == 0) {
            continue;
        }

        const i128 sp = prefix.prime_sum(R) - prefix.prime_sum(L - 1);
        const i128 group = tri(q) * (sp - 2 * static_cast<i128>(cnt));
        add_mod_i128(group);
    }

    return answer;
}

u64 brute_prefix_formula(int N) {
    i128 s = 0;
    for (int n = 1; n <= N; ++n) {
        s += t_formula_small(n);
    }
    s %= static_cast<i128>(kMod);
    if (s < 0) {
        s += static_cast<i128>(kMod);
    }
    return static_cast<u64>(s);
}

void run_validations() {
    assert(t_bruteforce_graph(45) == 52);
    for (int n = 1; n <= 200; ++n) {
        assert(t_formula_small(n) == t_bruteforce_graph(n));
    }

    assert(solve(10) == 26);
    assert(solve(100) == 5282);
    assert(solve(1000) == brute_prefix_formula(1000));
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kN = 1'000'000'000'000ULL;
    std::cout << solve(kN) << '\n';
    return 0;
}
