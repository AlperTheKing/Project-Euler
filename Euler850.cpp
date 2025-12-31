#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

using namespace std;

namespace {

constexpr size_t kMaxK = 24;
constexpr int kMaxExpLimit = 64;

struct Context {
    uint64_t N = 0;
    uint64_t mod2 = 0;
    int max_exp = 0;
    size_t k_count = 0;
    const vector<uint32_t>* primes = nullptr;
    const vector<int8_t>* mu = nullptr;
    const vector<vector<int>>* a_table = nullptr;
};

static uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t mod) {
    return static_cast<uint64_t>((__int128)a * b % mod);
}

static uint64_t isqrt_u64(uint64_t x) {
    long double approx = sqrtl(static_cast<long double>(x));
    uint64_t r = static_cast<uint64_t>(approx);
    while ((r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

static vector<int8_t> mobius_sieve(int max_n, vector<uint32_t>& primes_out) {
    vector<int> lp(max_n + 1, 0);
    vector<int8_t> mu(max_n + 1, 0);
    mu[1] = 1;
    for (int i = 2; i <= max_n; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes_out.push_back(static_cast<uint32_t>(i));
            mu[i] = -1;
        }
        for (int p : primes_out) {
            int64_t v = static_cast<int64_t>(i) * p;
            if (v > max_n) break;
            lp[v] = p;
            if (i % p == 0) {
                mu[v] = 0;
                break;
            }
            mu[v] = static_cast<int8_t>(-mu[i]);
        }
    }
    return mu;
}

// Counts squarefree s <= x with gcd(s, r)=1 where r is squarefree and encoded by its divisors/signs.
static uint64_t count_squarefree_coprime(uint64_t x,
                                         const vector<uint32_t>& primes_in_r,
                                         const vector<uint64_t>& divs,
                                         const vector<int8_t>& signs,
                                         const vector<int8_t>& mu) {
    if (x == 0) return 0;
    uint64_t limit = isqrt_u64(x);
    int64_t res = 0;
    for (uint64_t d = 1; d <= limit; ++d) {
        int8_t mu_d = mu[d];
        if (mu_d == 0) continue;
        bool ok = true;
        for (uint32_t p : primes_in_r) {
            if (d % p == 0) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;
        uint64_t y = x / (d * d);
        int64_t phi = 0;
        for (size_t i = 0; i < divs.size(); ++i) {
            phi += signs[i] * static_cast<int64_t>(y / divs[i]);
        }
        res += mu_d * phi;
    }
    return static_cast<uint64_t>(res);
}

template <size_t KMAX>
static void dfs(int start_idx,
                uint64_t curr_u,
                vector<uint32_t>& primes_in_r,
                vector<uint64_t>& divs,
                vector<int8_t>& signs,
                const array<uint64_t, KMAX>& fvals,
                const Context& ctx,
                vector<uint64_t>& acc) {
    uint64_t x = ctx.N / curr_u;
    uint64_t q = count_squarefree_coprime(x, primes_in_r, divs, signs, *ctx.mu);
    uint64_t q_mod = q % ctx.mod2;
    for (size_t i = 0; i < ctx.k_count; ++i) {
        acc[i] = (acc[i] + mul_mod(fvals[i], q_mod, ctx.mod2)) % ctx.mod2;
    }

    const auto& primes = *ctx.primes;
    const auto& a_table = *ctx.a_table;
    for (int i = start_idx; i < static_cast<int>(primes.size()); ++i) {
        uint64_t p = primes[i];
        if (curr_u > ctx.N / (p * p)) break;
        size_t primes_size = primes_in_r.size();
        size_t divs_size = divs.size();
        primes_in_r.push_back(static_cast<uint32_t>(p));
        for (size_t j = 0; j < divs_size; ++j) {
            divs.push_back(divs[j] * p);
            signs.push_back(static_cast<int8_t>(-signs[j]));
        }

        array<uint64_t, kMaxExpLimit> p_pow_mod{};
        p_pow_mod[0] = 1 % ctx.mod2;
        for (int a = 1; a <= ctx.max_exp; ++a) {
            p_pow_mod[a] = mul_mod(p_pow_mod[a - 1], p, ctx.mod2);
        }

        uint64_t p_pow = p * p;
        int e = 2;
        while (curr_u <= ctx.N / p_pow) {
            uint64_t new_u = curr_u * p_pow;
            array<uint64_t, KMAX> new_fvals = fvals;
            for (size_t k_idx = 0; k_idx < ctx.k_count; ++k_idx) {
                int a = a_table[k_idx][e];
                new_fvals[k_idx] = mul_mod(new_fvals[k_idx], p_pow_mod[a], ctx.mod2);
            }
            dfs(i + 1, new_u, primes_in_r, divs, signs, new_fvals, ctx, acc);
            if (p_pow > ctx.N / p) break;
            p_pow *= p;
            ++e;
        }

        primes_in_r.resize(primes_size);
        divs.resize(divs_size);
        signs.resize(divs_size);
    }
}

template <size_t KMAX>
static void worker(const Context& ctx,
                   int prime_limit,
                   atomic<int>& next_idx,
                   vector<uint64_t>& acc) {
    const auto& primes = *ctx.primes;
    const auto& a_table = *ctx.a_table;
    while (true) {
        int idx = next_idx.fetch_add(1);
        if (idx >= prime_limit) break;
        uint64_t p = primes[idx];
        if (p * p > ctx.N) continue;

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

        array<uint64_t, KMAX> base_fvals{};
        for (size_t i = 0; i < ctx.k_count; ++i) base_fvals[i] = 1;

        array<uint64_t, kMaxExpLimit> p_pow_mod{};
        p_pow_mod[0] = 1 % ctx.mod2;
        for (int a = 1; a <= ctx.max_exp; ++a) {
            p_pow_mod[a] = mul_mod(p_pow_mod[a - 1], p, ctx.mod2);
        }

        uint64_t p_pow = p * p;
        int e = 2;
        while (p_pow <= ctx.N) {
            array<uint64_t, KMAX> new_fvals = base_fvals;
            for (size_t k_idx = 0; k_idx < ctx.k_count; ++k_idx) {
                int a = a_table[k_idx][e];
                new_fvals[k_idx] = mul_mod(new_fvals[k_idx], p_pow_mod[a], ctx.mod2);
            }
            dfs(idx + 1, p_pow, primes_in_r, divs, signs, new_fvals, ctx, acc);
            if (p_pow > ctx.N / p) break;
            p_pow *= p;
            ++e;
        }
    }
}

static int max_power_of_two(uint64_t n) {
    int e = 0;
    while (n >= 2) {
        n >>= 1;
        ++e;
    }
    return e;
}

// Computes T(N)=2*S(N) mod 2*mod by summing over powerful parts and squarefree coprime counts.
static uint64_t compute_T_mod(uint64_t N, uint64_t mod, bool use_threads) {
    uint64_t mod2 = mod * 2;
    int max_e = max_power_of_two(N);
    assert(max_e < kMaxExpLimit);
    int large_start = (max_e % 2 == 1) ? max_e : max_e + 1;
    vector<int> small_ks;
    for (int k = 3; k < large_start; k += 2) small_ks.push_back(k);
    size_t k_small = small_ks.size();
    size_t k_count = k_small + 1;
    assert(k_count <= kMaxK);

    vector<vector<int>> a_table(k_count, vector<int>(max_e + 1, 0));
    for (size_t idx = 0; idx < k_small; ++idx) {
        int k = small_ks[idx];
        for (int e = 2; e <= max_e; ++e) {
            int c = (e + k - 1) / k;
            a_table[idx][e] = e - c;
        }
    }
    for (int e = 2; e <= max_e; ++e) {
        a_table[k_small][e] = e - 1;
    }

    uint64_t sqrt_n = isqrt_u64(N);
    vector<uint32_t> primes;
    vector<int8_t> mu = mobius_sieve(static_cast<int>(sqrt_n), primes);

    Context ctx;
    ctx.N = N;
    ctx.mod2 = mod2;
    ctx.max_exp = max_e;
    ctx.k_count = k_count;
    ctx.primes = &primes;
    ctx.mu = &mu;
    ctx.a_table = &a_table;

    vector<uint64_t> acc(k_count, 0);

    vector<uint32_t> primes_in_r;
    vector<uint64_t> divs = {1};
    vector<int8_t> signs = {1};
    array<uint64_t, kMaxK> base_fvals{};
    for (size_t i = 0; i < k_count; ++i) base_fvals[i] = 1;
    uint64_t q_root = count_squarefree_coprime(N, primes_in_r, divs, signs, mu);
    uint64_t q_root_mod = q_root % mod2;
    for (size_t i = 0; i < k_count; ++i) {
        acc[i] = (acc[i] + q_root_mod) % mod2;
    }

    int prime_limit = 0;
    while (prime_limit < static_cast<int>(primes.size()) &&
           static_cast<uint64_t>(primes[prime_limit]) * primes[prime_limit] <= N) {
        ++prime_limit;
    }

    int hw_threads = static_cast<int>(thread::hardware_concurrency());
    if (!use_threads || hw_threads <= 1 || prime_limit == 0) {
        atomic<int> next_idx(0);
        worker<kMaxK>(ctx, prime_limit, next_idx, acc);
    } else {
        int thread_count = hw_threads;
        atomic<int> next_idx(0);
        vector<thread> threads;
        vector<vector<uint64_t>> partials(thread_count, vector<uint64_t>(k_count, 0));
        threads.reserve(thread_count);
        for (int t = 0; t < thread_count; ++t) {
            threads.emplace_back([&ctx, prime_limit, &next_idx, &partials, t]() {
                worker<kMaxK>(ctx, prime_limit, next_idx, partials[t]);
            });
        }
        for (auto& th : threads) th.join();
        for (int t = 0; t < thread_count; ++t) {
            for (size_t i = 0; i < k_count; ++i) {
                acc[i] = (acc[i] + partials[t][i]) % mod2;
            }
        }
    }

    uint64_t h_mod = N % mod2;
    for (size_t i = 0; i < k_small; ++i) {
        h_mod = (h_mod + acc[i]) % mod2;
    }
    uint64_t f_large = acc[k_small];
    uint64_t large_start_k = static_cast<uint64_t>(max(large_start, 3));
    uint64_t count_large = 0;
    if (large_start_k <= N) {
        count_large = (N - large_start_k) / 2 + 1;
    }
    h_mod = (h_mod + mul_mod(count_large % mod2, f_large, mod2)) % mod2;

    __int128 sum_n = static_cast<__int128>(N) * (N + 1) / 2;
    uint64_t sum_n_mod = static_cast<uint64_t>(sum_n % mod2);
    uint64_t odd_count = (N + 1) / 2;
    uint64_t odd_mod = odd_count % mod2;
    uint64_t t_mod = mul_mod(odd_mod, sum_n_mod, mod2);
    if (t_mod >= h_mod) t_mod -= h_mod;
    else t_mod = t_mod + mod2 - h_mod;
    return t_mod;
}

static uint64_t pow_int(uint64_t base, int exp) {
    uint64_t res = 1;
    for (int i = 0; i < exp; ++i) res *= base;
    return res;
}

static uint64_t brute_T_small(int N) {
    vector<int> spf(N + 1, 0);
    for (int i = 2; i <= N; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            if (static_cast<int64_t>(i) * i <= N) {
                for (int j = i * i; j <= N; j += i) {
                    if (spf[j] == 0) spf[j] = i;
                }
            }
        }
    }
    vector<vector<pair<int, int>>> factors(N + 1);
    for (int n = 2; n <= N; ++n) {
        int x = n;
        while (x > 1) {
            int p = spf[x];
            int e = 0;
            while (x % p == 0) {
                x /= p;
                ++e;
            }
            factors[n].push_back({p, e});
        }
    }
    uint64_t T = 0;
    for (int k = 1; k <= N; k += 2) {
        for (int n = 1; n <= N; ++n) {
            uint64_t m = 1;
            for (const auto& pe : factors[n]) {
                int c = (pe.second + k - 1) / k;
                m *= pow_int(pe.first, c);
            }
            T += static_cast<uint64_t>(n - n / m);
        }
    }
    return T;
}

static uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1 % mod;
    uint64_t x = base % mod;
    uint64_t e = exp;
    while (e > 0) {
        if (e & 1) res = mul_mod(res, x, mod);
        x = mul_mod(x, x, mod);
        e >>= 1;
    }
    return res;
}

static void validate_formula() {
    for (int n = 2; n <= 20; ++n) {
        vector<pair<int, int>> fac;
        int x = n;
        for (int p = 2; p * p <= x; ++p) {
            if (x % p != 0) continue;
            int e = 0;
            while (x % p == 0) {
                x /= p;
                ++e;
            }
            fac.push_back({p, e});
        }
        if (x > 1) fac.push_back({x, 1});
        for (int k : {1, 3, 5, 7, 9}) {
            uint64_t m = 1;
            for (const auto& pe : fac) {
                int c = (pe.second + k - 1) / k;
                m *= pow_int(pe.first, c);
            }
            uint64_t expected2 = static_cast<uint64_t>(n - n / m);
            uint64_t sum_rem = 0;
            for (int i = 1; i <= n; ++i) {
                sum_rem += mod_pow(i, k, n);
            }
            uint64_t lhs = 2 * sum_rem;
            uint64_t rhs = static_cast<uint64_t>(n) * expected2;
            if (lhs != rhs) {
                cerr << "Validation failed for n=" << n << " k=" << k << '\n';
                exit(1);
            }
        }
    }
}

static void validate_small() {
    const uint64_t expected_10 = 201;
    const uint64_t expected_1000 = 247375608;
    uint64_t t10 = brute_T_small(10);
    uint64_t t1000 = brute_T_small(1000);
    if (t10 != expected_10) {
        cerr << "Validation failed for N=10: got " << t10 << '\n';
        exit(1);
    }
    if (t1000 != expected_1000) {
        cerr << "Validation failed for N=1000: got " << t1000 << '\n';
        exit(1);
    }
    const uint64_t mod = 977676779ULL;
    uint64_t t10_mod = compute_T_mod(10, mod, false);
    uint64_t t1000_mod = compute_T_mod(1000, mod, false);
    if (t10_mod != expected_10) {
        cerr << "Algorithm check failed for N=10: got " << t10_mod << '\n';
        exit(1);
    }
    if (t1000_mod != expected_1000) {
        cerr << "Algorithm check failed for N=1000: got " << t1000_mod << '\n';
        exit(1);
    }
}

}  // namespace

int main() {
    validate_formula();
    validate_small();

    const uint64_t N = 33557799775533ULL;
    const uint64_t MOD = 977676779ULL;

    uint64_t t_mod = compute_T_mod(N, MOD, true);
    uint64_t answer = (t_mod - (t_mod & 1)) / 2;
    cout << answer % MOD << '\n';
    return 0;
}
