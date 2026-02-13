#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct PrimeData {
    int exp = 0;
    std::vector<u64> phi;
    std::vector<u64> ord;
};

std::vector<int> sieve_primes(int limit) {
    std::vector<bool> is_prime(limit + 1, true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; 1LL * p * p <= limit; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[q] = false;
        }
    }
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

std::vector<std::pair<u64, int>> factorize_u64(u64 x, const std::vector<int>& primes) {
    std::vector<std::pair<u64, int>> out;
    u64 n = x;
    for (int p : primes) {
        if (1ULL * p * p > n) {
            break;
        }
        if (n % static_cast<u64>(p) != 0) {
            continue;
        }
        int e = 0;
        while (n % static_cast<u64>(p) == 0) {
            n /= static_cast<u64>(p);
            ++e;
        }
        out.push_back({static_cast<u64>(p), e});
    }
    if (n > 1) {
        out.push_back({n, 1});
    }
    return out;
}

void add_factorization(std::map<u64, int>& mp, u64 x, const std::vector<int>& primes) {
    if (x <= 1) {
        return;
    }
    const auto f = factorize_u64(x, primes);
    for (const auto& [p, e] : f) {
        mp[p] += e;
    }
}

u64 mod_pow(u64 a, u64 e, u64 mod) {
    if (mod == 1) {
        return 0;
    }
    u64 base = a % mod;
    u64 exp = e;
    u64 res = 1 % mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            res = static_cast<u64>((u128)res * base % mod);
        }
        base = static_cast<u64>((u128)base * base % mod);
        exp >>= 1ULL;
    }
    return res;
}

u64 ipow_u64(u64 a, int e) {
    u64 r = 1;
    for (int i = 0; i < e; ++i) {
        r *= a;
    }
    return r;
}

u64 lcm_u64(u64 a, u64 b) {
    if (a == 0 || b == 0) {
        return 0;
    }
    const u64 g = std::gcd(a, b);
    return static_cast<u64>((u128)(a / g) * b);
}

const std::vector<std::pair<u64, int>>& factor_cached(
    u64 x,
    const std::vector<int>& primes,
    std::unordered_map<u64, std::vector<std::pair<u64, int>>>& cache) {
    auto it = cache.find(x);
    if (it != cache.end()) {
        return it->second;
    }
    auto [jt, _] = cache.emplace(x, factorize_u64(x, primes));
    return jt->second;
}

u64 order_prime_power(
    u64 base,
    u64 p,
    int k,
    const std::vector<int>& primes,
    std::unordered_map<u64, std::vector<std::pair<u64, int>>>& factor_cache) {
    const u64 mod = ipow_u64(p, k);
    const u64 phi = ipow_u64(p, k - 1) * (p - 1);

    std::vector<std::pair<u64, int>> fac = factor_cached(p - 1, primes, factor_cache);
    if (k > 1) {
        bool found = false;
        for (auto& [q, e] : fac) {
            if (q == p) {
                e += k - 1;
                found = true;
                break;
            }
        }
        if (!found) {
            fac.push_back({p, k - 1});
        }
    }

    u64 ord = phi;
    for (const auto& [q, e] : fac) {
        for (int i = 0; i < e; ++i) {
            if (ord % q != 0) {
                break;
            }
            const u64 cand = ord / q;
            if (mod_pow(base, cand, mod) == 1) {
                ord = cand;
            } else {
                break;
            }
        }
    }
    return ord;
}

u64 swaps_fourth_power_dims(
    int n,
    int m,
    const std::vector<int>& primes,
    std::unordered_map<u64, std::vector<std::pair<u64, int>>>& factor_cache,
    std::unordered_map<int, std::vector<std::pair<u64, int>>>& t_factor_cache) {
    const u64 uu = static_cast<u64>(n) * static_cast<u64>(m);
    const u64 N = ipow_u64(uu, 4);
    const u64 M = N - 1;
    const u64 g = ipow_u64(static_cast<u64>(m), 4);

    std::vector<std::pair<u64, int>> factors;
    auto it_t = t_factor_cache.find(static_cast<int>(uu));
    if (it_t != t_factor_cache.end()) {
        factors = it_t->second;
    } else {
        std::map<u64, int> acc;
        add_factorization(acc, uu - 1, primes);
        add_factorization(acc, uu + 1, primes);
        add_factorization(acc, uu * uu + 1, primes);

        factors.reserve(acc.size());
        for (const auto& [p, e] : acc) {
            factors.push_back({p, e});
        }
        t_factor_cache.emplace(static_cast<int>(uu), factors);
    }

    std::vector<PrimeData> data;
    data.reserve(factors.size());

    for (const auto& [p, e] : factors) {
        PrimeData pd;
        pd.exp = e;
        pd.phi.assign(static_cast<std::size_t>(e + 1), 1);
        pd.ord.assign(static_cast<std::size_t>(e + 1), 1);

        for (int k = 1; k <= e; ++k) {
            pd.phi[static_cast<std::size_t>(k)] = ipow_u64(p, k - 1) * (p - 1);
            pd.ord[static_cast<std::size_t>(k)] =
                order_prime_power(g, p, k, primes, factor_cache);
        }

        data.push_back(std::move(pd));
    }

    u64 cycles_mod_set = 0;

    std::function<void(std::size_t, u64, u64)> dfs = [&](std::size_t idx, u64 phi_cur, u64 ord_cur) {
        if (idx == data.size()) {
            cycles_mod_set += phi_cur / ord_cur;
            return;
        }

        const PrimeData& pd = data[idx];
        for (int k = 0; k <= pd.exp; ++k) {
            const u64 phi_next = static_cast<u64>((u128)phi_cur * pd.phi[static_cast<std::size_t>(k)]);
            const u64 ord_next = lcm_u64(ord_cur, pd.ord[static_cast<std::size_t>(k)]);
            dfs(idx + 1, phi_next, ord_next);
        }
    };

    dfs(0, 1, 1);

    return N - 1 - cycles_mod_set;
}

u64 brute_swaps(u64 n, u64 m) {
    const u64 N = n * m;
    std::vector<char> vis(static_cast<std::size_t>(N), 0);

    auto next = [&](u64 x) -> u64 {
        if (x == N - 1) {
            return N - 1;
        }
        return static_cast<u64>((u128)x * m % (N - 1));
    };

    u64 cycles = 0;
    for (u64 i = 0; i < N; ++i) {
        if (vis[static_cast<std::size_t>(i)]) {
            continue;
        }
        ++cycles;
        u64 x = i;
        while (!vis[static_cast<std::size_t>(x)]) {
            vis[static_cast<std::size_t>(x)] = 1;
            x = next(x);
        }
    }

    return N - cycles;
}

void print_u128(u128 x) {
    if (x == 0) {
        std::cout << 0;
        return;
    }
    std::string s;
    u128 v = x;
    while (v > 0) {
        const int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }
    std::reverse(s.begin(), s.end());
    std::cout << s;
}

void validate(
    const std::vector<int>& primes,
    std::unordered_map<u64, std::vector<std::pair<u64, int>>>& factor_cache,
    std::unordered_map<int, std::vector<std::pair<u64, int>>>& t_factor_cache) {
    assert(brute_swaps(3, 4) == 8);

    u64 sum_small = 0;
    for (int n = 2; n <= 100; ++n) {
        for (int m = n; m <= 100; ++m) {
            sum_small += brute_swaps(static_cast<u64>(n), static_cast<u64>(m));
        }
    }
    assert(sum_small == 12'578'833ULL);

    for (int n = 2; n <= 5; ++n) {
        for (int m = n; m <= 5; ++m) {
            const u64 nn = ipow_u64(static_cast<u64>(n), 4);
            const u64 mm = ipow_u64(static_cast<u64>(m), 4);
            const u64 brute = brute_swaps(nn, mm);
            const u64 fast = swaps_fourth_power_dims(n, m, primes, factor_cache, t_factor_cache);
            assert(brute == fast);
        }
    }
}

}  // namespace

int main() {
    const std::vector<int> primes = sieve_primes(100'000);
    std::unordered_map<u64, std::vector<std::pair<u64, int>>> factor_cache;
    std::unordered_map<int, std::vector<std::pair<u64, int>>> t_factor_cache;

    validate(primes, factor_cache, t_factor_cache);

    u128 answer = 0;
    for (int n = 2; n <= 100; ++n) {
        for (int m = n; m <= 100; ++m) {
            answer += swaps_fourth_power_dims(n, m, primes, factor_cache, t_factor_cache);
        }
    }

    print_u128(answer);
    std::cout << '\n';
    return 0;
}
