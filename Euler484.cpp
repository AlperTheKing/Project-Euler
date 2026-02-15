#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using namespace std;

namespace {

using u128 = unsigned __int128;

struct Context {
    uint64_t N = 0;
    const vector<uint32_t>* primes = nullptr;
    const vector<int8_t>* mu = nullptr;
    const vector<uint32_t>* sqf_d = nullptr;
};

static uint64_t isqrt_u64(uint64_t x) {
    long double approx = sqrtl(static_cast<long double>(x));
    uint64_t r = static_cast<uint64_t>(approx);
    while ((r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

static vector<int8_t> mobius_sieve(int max_n, vector<uint32_t>& primes_out) {
    vector<int8_t> mu(max_n + 1, 0);
    vector<uint8_t> is_comp(max_n + 1, 0);
    mu[1] = 1;
    for (int i = 2; i <= max_n; ++i) {
        if (!is_comp[i]) {
            primes_out.push_back(static_cast<uint32_t>(i));
            mu[i] = -1;
        }
        for (uint32_t p : primes_out) {
            uint64_t v = static_cast<uint64_t>(i) * p;
            if (v > static_cast<uint64_t>(max_n)) break;
            is_comp[static_cast<size_t>(v)] = 1;
            if (i % static_cast<int>(p) == 0) {
                mu[static_cast<size_t>(v)] = 0;
                break;
            }
            mu[static_cast<size_t>(v)] = static_cast<int8_t>(-mu[i]);
        }
    }
    return mu;
}

// Counts squarefree s <= x with gcd(s, r)=1 where r is squarefree and encoded by its divisors/signs.
static uint64_t count_squarefree_coprime(uint64_t x,
                                         const vector<uint32_t>& primes_in_r,
                                         const vector<uint64_t>& divs,
                                         const vector<int8_t>& signs,
                                         const vector<int8_t>& mu,
                                         const vector<uint32_t>& sqf_d) {
    if (x == 0) return 0;
    uint64_t limit = isqrt_u64(x);
    auto end_it = upper_bound(sqf_d.begin(), sqf_d.end(), static_cast<uint32_t>(limit));
    const size_t prime_cnt = primes_in_r.size();

    int64_t res = 0;
    uint64_t last_y = numeric_limits<uint64_t>::max();
    int64_t last_phi = 0;

    for (auto it = sqf_d.begin(); it != end_it; ++it) {
        uint64_t d = *it;
        bool ok = true;
        if (prime_cnt == 1) {
            ok = (d % primes_in_r[0] != 0);
        } else if (prime_cnt == 2) {
            ok = (d % primes_in_r[0] != 0) && (d % primes_in_r[1] != 0);
        } else if (prime_cnt == 3) {
            ok = (d % primes_in_r[0] != 0) && (d % primes_in_r[1] != 0) &&
                 (d % primes_in_r[2] != 0);
        } else if (prime_cnt == 4) {
            ok = (d % primes_in_r[0] != 0) && (d % primes_in_r[1] != 0) &&
                 (d % primes_in_r[2] != 0) && (d % primes_in_r[3] != 0);
        } else {
            for (uint32_t p : primes_in_r) {
                if (d % p == 0) {
                    ok = false;
                    break;
                }
            }
        }
        if (!ok) continue;

        const int8_t mu_d = mu[static_cast<size_t>(d)];
        uint64_t y = x / (d * d);
        if (y != last_y) {
            int64_t phi = 0;
            for (size_t i = 0; i < divs.size(); ++i) {
                phi += signs[i] * static_cast<int64_t>(y / divs[i]);
            }
            last_y = y;
            last_phi = phi;
        }
        res += mu_d * last_phi;
    }
    return static_cast<uint64_t>(res);
}

// For prime power p^e (e>=2): g(p^e)=p^(e-1) unless p | e, then g=p^e.
static uint64_t g_prime_power(uint64_t p, int e, uint64_t p_pow) {
    if (static_cast<uint64_t>(e) % p == 0) return p_pow;
    return p_pow / p;
}

static void dfs(int start_idx,
                uint64_t curr_u,
                uint64_t curr_g,
                vector<uint32_t>& primes_in_r,
                vector<uint64_t>& divs,
                vector<int8_t>& signs,
                const Context& ctx,
                u128& acc) {
    // Each node corresponds to one powerful part; squarefree coprime factors are counted separately.
    uint64_t x = ctx.N / curr_u;
    uint64_t q = count_squarefree_coprime(x, primes_in_r, divs, signs, *ctx.mu, *ctx.sqf_d);
    acc += static_cast<u128>(curr_g) * q;

    const auto& primes = *ctx.primes;
    for (int i = start_idx; i < static_cast<int>(primes.size()); ++i) {
        uint64_t p = primes[i];
        uint64_t p2 = p * p;
        if (curr_u > ctx.N / p2) break;

        size_t primes_size = primes_in_r.size();
        size_t divs_size = divs.size();
        primes_in_r.push_back(static_cast<uint32_t>(p));
        for (size_t j = 0; j < divs_size; ++j) {
            divs.push_back(divs[j] * p);
            signs.push_back(static_cast<int8_t>(-signs[j]));
        }

        uint64_t p_pow = p2;
        int e = 2;
        while (curr_u <= ctx.N / p_pow) {
            uint64_t new_u = curr_u * p_pow;
            uint64_t g_p = g_prime_power(p, e, p_pow);
            uint64_t new_g = curr_g * g_p;
            dfs(i + 1, new_u, new_g, primes_in_r, divs, signs, ctx, acc);
            if (p_pow > ctx.N / p) break;
            p_pow *= p;
            ++e;
        }

        primes_in_r.resize(primes_size);
        divs.resize(divs_size);
        signs.resize(divs_size);
    }
}

static void worker(const Context& ctx, int prime_limit, atomic<int>& next_idx, u128& acc) {
    const auto& primes = *ctx.primes;
    while (true) {
        int idx = next_idx.fetch_add(1);
        if (idx >= prime_limit) break;
        uint64_t p = primes[idx];
        uint64_t p2 = p * p;
        if (p2 > ctx.N) continue;

        vector<uint32_t> primes_in_r;
        primes_in_r.reserve(8);
        primes_in_r.push_back(static_cast<uint32_t>(p));

        vector<uint64_t> divs;
        vector<int8_t> signs;
        divs.reserve(128);
        signs.reserve(128);
        divs.push_back(1);
        divs.push_back(p);
        signs.push_back(1);
        signs.push_back(-1);

        uint64_t p_pow = p2;
        int e = 2;
        while (p_pow <= ctx.N) {
            uint64_t g_p = g_prime_power(p, e, p_pow);
            dfs(idx + 1, p_pow, g_p, primes_in_r, divs, signs, ctx, acc);
            if (p_pow > ctx.N / p) break;
            p_pow *= p;
            ++e;
        }
    }
}

static u128 compute_sum(uint64_t N, bool use_threads) {
    uint64_t sqrt_n = isqrt_u64(N);
    vector<uint32_t> primes;
    vector<int8_t> mu = mobius_sieve(static_cast<int>(sqrt_n), primes);
    vector<uint32_t> sqf_d;
    sqf_d.reserve(static_cast<size_t>(sqrt_n * 7 / 10));
    for (uint32_t d = 1; d <= sqrt_n; ++d) {
        if (mu[static_cast<size_t>(d)] != 0) sqf_d.push_back(d);
    }

    Context ctx;
    ctx.N = N;
    ctx.primes = &primes;
    ctx.mu = &mu;
    ctx.sqf_d = &sqf_d;

    u128 total = 0;
    vector<uint32_t> primes_in_r;
    vector<uint64_t> divs = {1};
    vector<int8_t> signs = {1};
    uint64_t q_root = count_squarefree_coprime(N, primes_in_r, divs, signs, mu, sqf_d);
    total += q_root;

    int prime_limit = 0;
    while (prime_limit < static_cast<int>(primes.size())) {
        uint64_t p = primes[prime_limit];
        if (p * p > N) break;
        ++prime_limit;
    }

    int hw_threads = static_cast<int>(thread::hardware_concurrency());
    if (!use_threads || hw_threads <= 1 || prime_limit == 0) {
        atomic<int> next_idx(0);
        worker(ctx, prime_limit, next_idx, total);
    } else {
        int thread_count = hw_threads;
        atomic<int> next_idx(0);
        vector<thread> threads;
        vector<u128> partials(thread_count, 0);
        threads.reserve(thread_count);
        for (int t = 0; t < thread_count; ++t) {
            threads.emplace_back([&ctx, prime_limit, &next_idx, &partials, t]() {
                worker(ctx, prime_limit, next_idx, partials[t]);
            });
        }
        for (auto& th : threads) th.join();
        for (int t = 0; t < thread_count; ++t) {
            total += partials[t];
        }
    }
    return total;
}

static vector<uint32_t> generate_primes_upto(uint32_t limit) {
    vector<uint32_t> primes;
    vector<uint8_t> is_comp(static_cast<size_t>(limit + 1), 0);
    for (uint32_t i = 2; i <= limit; ++i) {
        if (!is_comp[static_cast<size_t>(i)]) {
            primes.push_back(i);
            if (static_cast<uint64_t>(i) * i <= limit) {
                for (uint64_t j = static_cast<uint64_t>(i) * i; j <= limit; j += i) {
                    is_comp[static_cast<size_t>(j)] = 1;
                }
            }
        }
    }
    return primes;
}

static u128 dfs_alt(int i0, uint64_t L0, const vector<uint32_t>& primes, const vector<uint64_t>& p2) {
    u128 res = 0;
    for (int i = i0; i < static_cast<int>(primes.size()); ++i) {
        const uint64_t q = p2[static_cast<size_t>(i)];
        uint64_t L = L0 / q;
        if (L == 0) break;

        const uint64_t p = primes[static_cast<size_t>(i)];
        uint64_t e = 1;
        uint64_t g = 1;
        while (L > 0) {
            const uint64_t gp = g;
            ++e;
            if (e != 1) {
                if (e == p) {
                    g *= q;
                    e = 0;
                } else {
                    g *= p;
                }
                const uint64_t c = g - gp;
                res += static_cast<u128>(c) * static_cast<u128>(L);
                if (L > q) {
                    res += static_cast<u128>(c) * dfs_alt(i + 1, L, primes, p2);
                }
            }
            L /= p;
        }
    }
    return res;
}

static u128 compute_sum_alt(uint64_t L) {
    const uint32_t limit = static_cast<uint32_t>(isqrt_u64(L));
    const vector<uint32_t> primes = generate_primes_upto(limit);
    vector<uint64_t> p2;
    p2.reserve(primes.size());
    for (uint32_t p : primes) p2.push_back(static_cast<uint64_t>(p) * p);
    return static_cast<u128>(L) + dfs_alt(0, L, primes, p2);
}

static vector<int> build_spf(int n) {
    vector<int> spf(n + 1, 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[i] != 0) continue;
        spf[i] = i;
        if (static_cast<int64_t>(i) * i > n) continue;
        for (int j = i * i; j <= n; j += i) {
            if (spf[j] == 0) spf[j] = i;
        }
    }
    return spf;
}

static uint64_t arithmetic_derivative(uint64_t n, const vector<int>& spf) {
    if (n == 1) return 0;
    uint64_t x = n;
    uint64_t sum = 0;
    while (x > 1) {
        int p = spf[static_cast<size_t>(x)];
        if (p == 0) p = static_cast<int>(x);
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        sum += static_cast<uint64_t>(e) * (n / p);
    }
    return sum;
}

static uint64_t gcd_formula_value(uint64_t n, const vector<int>& spf) {
    if (n == 1) return 1;
    uint64_t x = n;
    uint64_t g = 1;
    while (x > 1) {
        int p = spf[static_cast<size_t>(x)];
        if (p == 0) p = static_cast<int>(x);
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        uint64_t p_pow = 1;
        for (int i = 1; i < e; ++i) p_pow *= static_cast<uint64_t>(p);
        uint64_t term = p_pow;
        if (static_cast<uint64_t>(e) % static_cast<uint64_t>(p) == 0) {
            term *= static_cast<uint64_t>(p);
        }
        g *= term;
    }
    return g;
}

static uint64_t brute_sum(int n) {
    vector<int> spf = build_spf(n);
    uint64_t sum = 0;
    for (int i = 1; i <= n; ++i) {
        uint64_t deriv = arithmetic_derivative(static_cast<uint64_t>(i), spf);
        sum += gcd<uint64_t>(static_cast<uint64_t>(i), deriv);
    }
    return sum;
}

static void validate_formula() {
    const int limit = 5000;
    vector<int> spf = build_spf(limit);
    for (int n = 1; n <= limit; ++n) {
        uint64_t deriv = arithmetic_derivative(static_cast<uint64_t>(n), spf);
        uint64_t g = gcd<uint64_t>(static_cast<uint64_t>(n), deriv);
        uint64_t g_formula = gcd_formula_value(static_cast<uint64_t>(n), spf);
        if (g != g_formula) {
            cerr << "Formula check failed at n=" << n << '\n';
            exit(1);
        }
    }
}

static void validate_small() {
    const int n = 10000;
    uint64_t brute = brute_sum(n);
    u128 fast = compute_sum(static_cast<uint64_t>(n), false);
    if (fast != static_cast<u128>(brute)) {
        cerr << "Small-N check failed: expected " << brute << '\n';
        exit(1);
    }
    u128 alt = compute_sum_alt(static_cast<uint64_t>(n));
    if (alt != static_cast<u128>(brute)) {
        cerr << "Small-N alt check failed: expected " << brute
             << " got " << static_cast<uint64_t>(alt) << '\n';
        exit(1);
    }
}

static string to_string_u128(u128 value) {
    if (value == 0) return "0";
    string s;
    while (value > 0) {
        int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}

}  // namespace

int main() {
    validate_formula();
    validate_small();

    const uint64_t N = 5000000000000000ULL;
    u128 total = compute_sum_alt(N);
    if (total > 0) total -= 1;  // exclude n=1
    cout << to_string_u128(total) << '\n';
    return 0;
}
