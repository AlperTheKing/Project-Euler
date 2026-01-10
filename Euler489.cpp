#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <thread>
#include <utility>
#include <vector>

using u128 = unsigned __int128;

namespace {

struct PrimeSolutions {
    uint64_t mod = 1;
    std::vector<uint64_t> residues;
};

struct Precompute {
    std::vector<int> primes;
    std::vector<uint64_t> a_pow6;
    std::vector<uint64_t> b_term;
    std::vector<std::vector<std::pair<uint64_t, int>>> a_factors;
};

uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t mod) {
    return static_cast<uint64_t>((u128)a * b % mod);
}

uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1 % mod;
    uint64_t cur = base % mod;
    while (exp > 0) {
        if (exp & 1) result = mul_mod(result, cur, mod);
        cur = mul_mod(cur, cur, mod);
        exp >>= 1;
    }
    return result;
}

uint64_t mod_inv_prime(uint64_t a, uint64_t p) {
    return mod_pow(a, p - 2, p);
}

int64_t egcd(int64_t a, int64_t b, int64_t& x, int64_t& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int64_t x1 = 0;
    int64_t y1 = 0;
    int64_t g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return g;
}

uint64_t mod_inv_coprime(uint64_t a, uint64_t mod) {
    int64_t x = 0;
    int64_t y = 0;
    egcd(static_cast<int64_t>(a), static_cast<int64_t>(mod), x, y);
    int64_t res = x % static_cast<int64_t>(mod);
    if (res < 0) res += static_cast<int64_t>(mod);
    return static_cast<uint64_t>(res);
}

bool mod_sqrt(uint64_t n, uint64_t p, uint64_t& out) {
    if (p == 2) {
        out = n % p;
        return true;
    }
    n %= p;
    if (n == 0) {
        out = 0;
        return true;
    }
    if (mod_pow(n, (p - 1) / 2, p) != 1) return false;
    if (p % 4 == 3) {
        out = mod_pow(n, (p + 1) / 4, p);
        return true;
    }

    uint64_t q = p - 1;
    int s = 0;
    while ((q & 1) == 0) {
        q >>= 1;
        ++s;
    }

    uint64_t z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) ++z;

    uint64_t c = mod_pow(z, q, p);
    uint64_t x = mod_pow(n, (q + 1) / 2, p);
    uint64_t t = mod_pow(n, q, p);
    int m = s;
    while (t != 1) {
        int i = 1;
        uint64_t t2 = mul_mod(t, t, p);
        while (t2 != 1) {
            t2 = mul_mod(t2, t2, p);
            ++i;
        }
        uint64_t b = mod_pow(c, 1ULL << (m - i - 1), p);
        x = mul_mod(x, b, p);
        uint64_t b2 = mul_mod(b, b, p);
        t = mul_mod(t, b2, p);
        c = b2;
        m = i;
    }
    out = x;
    return true;
}

bool is_root_mod(uint64_t n, int a, int b, uint64_t mod) {
    u128 nn = n;
    u128 na = nn + static_cast<u128>(a);
    u128 val1 = nn * nn * nn + static_cast<u128>(b);
    u128 val2 = na * na * na + static_cast<u128>(b);
    return (val1 % mod) == 0 && (val2 % mod) == 0;
}

std::vector<int> sieve_primes(int limit) {
    std::vector<bool> is_prime(limit + 1, true);
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (!is_prime[i]) continue;
        primes.push_back(i);
        if (static_cast<int64_t>(i) * i > limit) continue;
        for (int j = i * i; j <= limit; j += i) {
            is_prime[j] = false;
        }
    }
    return primes;
}

std::vector<std::pair<uint64_t, int>> factorize(uint64_t n, const std::vector<int>& primes) {
    std::vector<std::pair<uint64_t, int>> out;
    uint64_t tmp = n;
    for (int p : primes) {
        if (static_cast<uint64_t>(p) * p > tmp) break;
        if (tmp % p != 0) continue;
        int e = 0;
        while (tmp % p == 0) {
            tmp /= p;
            ++e;
        }
        out.push_back({static_cast<uint64_t>(p), e});
    }
    if (tmp > 1) out.push_back({tmp, 1});
    return out;
}

void add_factor(std::vector<std::pair<uint64_t, int>>& factors, uint64_t p, int e) {
    for (auto& kv : factors) {
        if (kv.first == p) {
            kv.second += e;
            return;
        }
    }
    factors.push_back({p, e});
}

std::vector<uint64_t> solutions_mod_p(uint64_t p, int a, int b) {
    std::vector<uint64_t> sols;
    if (p <= 5 || (a % static_cast<int>(p)) == 0) {
        for (uint64_t n = 0; n < p; ++n) {
            if (is_root_mod(n, a, b, p)) sols.push_back(n);
        }
        return sols;
    }

    if (p == 3) return sols;

    uint64_t inv3 = mod_inv_prime(3, p);
    uint64_t inv2 = mod_inv_prime(2, p);
    uint64_t a_mod = static_cast<uint64_t>(a) % p;
    uint64_t D = (p + p - mul_mod(a_mod, a_mod, p)) % p;
    D = mul_mod(D, inv3, p);

    uint64_t sqrtD = 0;
    if (!mod_sqrt(D, p, sqrtD)) return sols;

    uint64_t roots[2] = {sqrtD, (p - sqrtD) % p};
    for (int i = 0; i < 2; ++i) {
        uint64_t s = roots[i];
        uint64_t tmp = (p + s - a_mod) % p;
        uint64_t n = mul_mod(tmp, inv2, p);
        if (is_root_mod(n, a, b, p)) {
            if (std::find(sols.begin(), sols.end(), n) == sols.end()) sols.push_back(n);
        }
    }
    return sols;
}

PrimeSolutions solve_prime_power(uint64_t p, int max_exp, int a, int b) {
    PrimeSolutions res;
    res.mod = 1;
    res.residues = {0};

    std::vector<uint64_t> sols = solutions_mod_p(p, a, b);
    if (sols.empty()) return res;

    uint64_t mod = p;
    for (int exp = 1; exp < max_exp; ++exp) {
        uint64_t next_mod = mod * p;
        std::vector<uint64_t> next;
        next.reserve(sols.size() * p);
        for (uint64_t r : sols) {
            for (uint64_t t = 0; t < p; ++t) {
                uint64_t n = r + t * mod;
                if (is_root_mod(n, a, b, next_mod)) {
                    next.push_back(n);
                }
            }
        }
        if (next.empty()) break;
        sols.swap(next);
        mod = next_mod;
    }

    res.mod = mod;
    res.residues = std::move(sols);
    return res;
}

uint64_t crt(uint64_t a, uint64_t mod_a, uint64_t b, uint64_t mod_b) {
    uint64_t inv = mod_inv_coprime(mod_a % mod_b, mod_b);
    uint64_t t = (b + mod_b - (a % mod_b)) % mod_b;
    uint64_t k = mul_mod(t, inv, mod_b);
    return static_cast<uint64_t>(a + static_cast<u128>(mod_a) * k);
}

uint64_t compute_G(int a, int b, const Precompute& pre) {
    uint64_t R = pre.a_pow6[a] + pre.b_term[b];
    std::vector<std::pair<uint64_t, int>> factors = factorize(R, pre.primes);
    for (const auto& kv : pre.a_factors[a]) {
        add_factor(factors, kv.first, 3 * kv.second);
    }

    std::vector<PrimeSolutions> blocks;
    blocks.reserve(factors.size());
    for (const auto& kv : factors) {
        PrimeSolutions sol = solve_prime_power(kv.first, kv.second, a, b);
        if (sol.mod > 1) blocks.push_back(std::move(sol));
    }
    if (blocks.empty()) return 0;

    std::sort(blocks.begin(), blocks.end(),
              [](const PrimeSolutions& lhs, const PrimeSolutions& rhs) {
                  return lhs.residues.size() < rhs.residues.size();
              });

    uint64_t mod = 1;
    std::vector<uint64_t> residues = {0};
    for (const auto& block : blocks) {
        std::vector<uint64_t> next;
        next.reserve(residues.size() * block.residues.size());
        for (uint64_t r : residues) {
            for (uint64_t s : block.residues) {
                next.push_back(crt(r, mod, s, block.mod));
            }
        }
        mod *= block.mod;
        residues.swap(next);
    }

    return *std::min_element(residues.begin(), residues.end());
}

uint64_t compute_H(int m, int n, const Precompute& pre, unsigned threads) {
    const int total = m * n;
    if (threads == 0) threads = 1;
    std::atomic<int> index(0);
    std::vector<std::thread> pool;
    std::vector<uint64_t> partial(threads, 0);

    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            while (true) {
                int i = index.fetch_add(1);
                if (i >= total) break;
                int a = i / n + 1;
                int b = i % n + 1;
                partial[t] += compute_G(a, b, pre);
            }
        });
    }
    for (auto& th : pool) th.join();

    uint64_t sum = 0;
    for (uint64_t v : partial) sum += v;
    return sum;
}

Precompute build_precompute(int max_a, int max_b) {
    Precompute pre;
    pre.primes = sieve_primes(200000);
    pre.a_pow6.assign(max_a + 1, 0);
    pre.a_factors.assign(max_a + 1, {});
    pre.b_term.assign(max_b + 1, 0);

    for (int a = 1; a <= max_a; ++a) {
        u128 v = 1;
        for (int i = 0; i < 6; ++i) v *= static_cast<u128>(a);
        pre.a_pow6[a] = static_cast<uint64_t>(v);
        pre.a_factors[a] = factorize(static_cast<uint64_t>(a), pre.primes);
    }

    for (int b = 1; b <= max_b; ++b) {
        pre.b_term[b] = 27ULL * static_cast<uint64_t>(b) * static_cast<uint64_t>(b);
    }

    return pre;
}

} // namespace

int main() {
    const int max_a = 18;
    const int max_b = 1900;
    Precompute pre = build_precompute(max_a, max_b);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;

    std::cout << "--- Validation Checkpoints ---\n";
    uint64_t g11 = compute_G(1, 1, pre);
    std::cout << "G(1, 1) = " << g11 << (g11 == 5 ? " [PASS]" : " [FAIL]") << "\n";

    uint64_t h55 = compute_H(5, 5, pre, threads);
    std::cout << "H(5, 5) = " << h55 << (h55 == 128878ULL ? " [PASS]" : " [FAIL]") << "\n";

    uint64_t h1010 = compute_H(10, 10, pre, threads);
    std::cout << "H(10, 10) = " << h1010 << (h1010 == 32936544ULL ? " [PASS]" : " [FAIL]") << "\n";

    std::cout << "\n--- Final Solution ---\n";
    uint64_t result = compute_H(18, 1900, pre, threads);
    std::cout << "H(18, 1900) = " << result << "\n";
    return 0;
}
