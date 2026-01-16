#include <algorithm>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using int64 = long long;

constexpr int64 kMod = 1000000007LL;
constexpr int64 kMainN = 1000000000000LL;
constexpr int kMainFact = 200;
constexpr int kThreshold = 10000000;  // Precompute sum_{n<=x} tau(n) up to 1e7.

int64 mod_norm(int64 v) {
    v %= kMod;
    if (v < 0) v += kMod;
    return v;
}

int64 mod_mul(int64 a, int64 b) {
    return static_cast<int64>((static_cast<__int128>(a) * b) % kMod);
}

struct SplitMix64 {
    std::size_t operator()(std::uint64_t x) const noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<std::size_t>(x ^ (x >> 31));
    }
};

std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int i = 2; i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

int64 factorial_prime_exp(int n, int p) {
    int64 exp = 0;
    while (n > 0) {
        n /= p;
        exp += n;
    }
    return exp;
}

int64 mod_pow(int64 base, int64 exp) {
    int64 res = 1 % kMod;
    int64 cur = base % kMod;
    while (exp > 0) {
        if (exp & 1LL) res = mod_mul(res, cur);
        cur = mod_mul(cur, cur);
        exp >>= 1LL;
    }
    return res;
}

struct PrimeCoeff {
    int p = 0;
    int64 g1 = 0;
};

struct CoeffData {
    std::vector<PrimeCoeff> coeffs;
    int64 constant = 1;
};

// Normalize f(k) by C = f(1) so F(k)=f(k)/C has F(p^b)=1+alpha*b and F = tau * g with squarefree g.
CoeffData build_coeffs(int fact_n) {
    std::vector<int> primes = primes_up_to(fact_n);
    CoeffData data;
    data.coeffs.reserve(primes.size());
    for (int p : primes) {
        int64 a = factorial_prime_exp(fact_n, p);
        int64 A = a + 1;
        int64 B = A * (A + 1) / 2;
        int64 B_mod = mod_norm(B);
        data.constant = mod_mul(data.constant, B_mod);

        int64 alpha = mod_mul(A % kMod, mod_pow(B_mod, kMod - 2));
        int64 g1 = mod_norm(alpha - 1);
        data.coeffs.push_back({p, g1});
    }
    return data;
}

struct GPrefixCalculator {
    std::vector<PrimeCoeff> coeffs;
    std::vector<std::unordered_map<int64, int64, SplitMix64>> memo;

    explicit GPrefixCalculator(std::vector<PrimeCoeff> coeffs_in)
        : coeffs(std::move(coeffs_in)), memo(coeffs.size()) {}

    int64 prefix(int64 limit) {
        return prefix_rec(0, limit);
    }

    int64 prefix_rec(std::size_t idx, int64 limit) {
        if (limit <= 0) return 0;
        if (idx == coeffs.size()) return 1;
        auto& cache = memo[idx];
        auto it = cache.find(limit);
        if (it != cache.end()) return it->second;

        const PrimeCoeff& coeff = coeffs[idx];
        int64 res = prefix_rec(idx + 1, limit);
        if (limit >= coeff.p) {
            res += mod_mul(coeff.g1, prefix_rec(idx + 1, limit / coeff.p));
            res %= kMod;
        }
        cache.emplace(limit, res);
        return res;
    }
};

std::vector<int> build_tau_prefix(int limit) {
    std::vector<int> tau(limit + 1, 0);
    for (int i = 1; i <= limit; ++i) {
        for (int j = i; j <= limit; j += i) {
            ++tau[j];
        }
    }
    std::vector<int> prefix(limit + 1, 0);
    int64 sum = 0;
    for (int i = 1; i <= limit; ++i) {
        sum += tau[i];
        if (sum >= kMod) sum -= kMod;
        prefix[i] = static_cast<int>(sum);
    }
    return prefix;
}

int64 divisor_summatory_mod(int64 n) {
    int64 res = 0;
    int64 i = 1;
    while (i <= n) {
        int64 q = n / i;
        int64 r = n / q;
        int64 cnt = r - i + 1;
        int64 term = mod_mul(q % kMod, cnt % kMod);
        res += term;
        if (res >= kMod) res -= kMod;
        i = r + 1;
    }
    return res;
}

std::vector<int64> collect_large_values(int64 n, int threshold) {
    std::vector<int64> values;
    int64 l = 1;
    while (l <= n) {
        int64 v = n / l;
        int64 r = n / v;
        if (v <= threshold) break;
        values.push_back(v);
        l = r + 1;
    }
    return values;
}

std::vector<int64> compute_large_T(const std::vector<int64>& values) {
    std::vector<int64> result(values.size(), 0);
    if (values.empty()) return result;

    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    if (values.size() < threads) threads = static_cast<unsigned int>(values.size());

    std::vector<std::thread> workers;
    workers.reserve(threads);
    const std::size_t chunk = (values.size() + threads - 1) / threads;

    for (unsigned int t = 0; t < threads; ++t) {
        std::size_t start = t * chunk;
        std::size_t end = std::min(values.size(), start + chunk);
        if (start >= end) break;
        workers.emplace_back([&, start, end]() {
            for (std::size_t i = start; i < end; ++i) {
                result[i] = divisor_summatory_mod(values[i]);
            }
        });
    }

    for (auto& th : workers) th.join();
    return result;
}

int64 compute_D(int fact_n, int64 n, const std::vector<int>& tau_prefix, int threshold) {
    std::vector<int64> large_values;
    std::vector<int64> large_T;
    if (n > threshold) {
        large_values = collect_large_values(n, threshold);
        large_T = compute_large_T(large_values);
    }

    CoeffData coeff_data = build_coeffs(fact_n);
    GPrefixCalculator gcalc(std::move(coeff_data.coeffs));

    int64 sum_f = 0;
    int64 l = 1;
    std::size_t large_idx = 0;
    while (l <= n) {
        int64 v = n / l;
        int64 r = n / v;

        int64 sum_g = gcalc.prefix(r) - gcalc.prefix(l - 1);
        if (sum_g < 0) sum_g += kMod;

        int64 t_val = 0;
        if (v <= threshold) {
            t_val = tau_prefix[static_cast<std::size_t>(v)];
        } else {
            t_val = large_T[large_idx++];
        }

        sum_f += mod_mul(sum_g, t_val);
        sum_f %= kMod;

        l = r + 1;
    }

    return mod_mul(coeff_data.constant, sum_f);
}

}  // namespace

int main() {
    const std::vector<int> tau_prefix = build_tau_prefix(kThreshold);

    struct Check {
        int fact_n;
        int64 n;
        int64 expected;
    };

    const Check checks[] = {
        {3, 100, 3398},
        {4, 1000000LL, 268882292},
    };

    for (const auto& chk : checks) {
        int64 got = compute_D(chk.fact_n, chk.n, tau_prefix, kThreshold);
        if (got != mod_norm(chk.expected)) {
            std::cerr << "Validation failure: D(" << chk.fact_n << "!, " << chk.n
                      << ") mod 1e9+7 = " << got
                      << ", expected " << mod_norm(chk.expected) << '\n';
            return 1;
        }
    }

    int64 answer = compute_D(kMainFact, kMainN, tau_prefix, kThreshold);
    std::cout << answer << '\n';
    return 0;
}
