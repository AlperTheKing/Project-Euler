#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <utility>
#include <unistd.h>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 MOD = 999'999'893ULL;

u64 add_mod(u64 a, u64 b) {
    a += b;
    if (a >= MOD) {
        a -= MOD;
    }
    return a;
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_u64(u64 base, int exp) {
    u64 result = 1;
    while (exp > 0) {
        if (exp & 1) {
            result *= base;
        }
        base *= base;
        exp >>= 1;
    }
    return result;
}

std::pair<std::vector<int>, std::vector<int>> build_spf_and_primes(int limit) {
    std::vector<int> spf(limit + 1, 0);
    std::vector<int> primes;
    primes.reserve(limit / 10);

    for (int i = 2; i <= limit; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const long long v = static_cast<long long>(p) * i;
            if (v > limit || p > spf[i]) {
                break;
            }
            spf[static_cast<int>(v)] = p;
        }
    }
    return {std::move(spf), std::move(primes)};
}

std::vector<std::pair<u32, int>> factorize_u32(u32 n, const std::vector<int>& spf) {
    std::vector<std::pair<u32, int>> factors;
    while (n > 1) {
        const u32 p = static_cast<u32>(spf[n]);
        int e = 0;
        do {
            n /= p;
            ++e;
        } while (n > 1 && static_cast<u32>(spf[n]) == p);
        factors.push_back({p, e});
    }
    return factors;
}

std::vector<u32> divisors_from_factorization(const std::vector<std::pair<u32, int>>& factors, bool sort_result) {
    std::vector<u32> divisors{1};
    for (const auto& [p, e] : factors) {
        const std::size_t current = divisors.size();
        u32 pe = 1;
        for (int i = 1; i <= e; ++i) {
            pe *= p;
            for (std::size_t j = 0; j < current; ++j) {
                divisors.push_back(divisors[j] * pe);
            }
        }
    }
    if (sort_result) {
        std::sort(divisors.begin(), divisors.end());
    }
    return divisors;
}

std::pair<u64, u64> fib_pair_mod(u64 n, u64 mod) {
    if (n == 0) {
        return {0ULL, 1ULL % mod};
    }
    const auto [a, b] = fib_pair_mod(n >> 1, mod);

    const u64 two_b = (2ULL * b) % mod;
    const u64 two_b_minus_a = (two_b + mod - a) % mod;
    const u64 c = mul_mod(a, two_b_minus_a, mod);
    const u64 d = (mul_mod(a, a, mod) + mul_mod(b, b, mod)) % mod;

    if (n & 1ULL) {
        return {d, (c + d) % mod};
    }
    return {c, d};
}

bool is_identity_power_mod_prime(u32 n, u32 p) {
    const auto [fn, fn1] = fib_pair_mod(n, p);
    return fn == 0ULL && fn1 == 1ULL;
}

u32 period_prime(u32 p, const std::vector<int>& spf) {
    if (p == 2U) {
        return 3U;
    }
    if (p == 5U) {
        return 20U;
    }

    const int legendre5 = (p % 5U == 1U || p % 5U == 4U) ? 1 : -1;
    u32 order = (legendre5 == 1) ? (p - 1U) : (2U * (p + 1U));

    const auto fac = factorize_u32(order, spf);
    for (const auto& [q, _] : fac) {
        while (order % q == 0U && is_identity_power_mod_prime(order / q, p)) {
            order /= q;
        }
    }
    return order;
}

int vp_limited(u64 x, u32 p, int cap) {
    int v = 0;
    while (v < cap && x % p == 0ULL) {
        x /= p;
        ++v;
    }
    return v;
}

u64 count_fixed_prime_power(u32 p, int e, u32 d) {
    const int cap = 2 * e;
    const u64 mod = pow_u64(p, cap);

    const auto [fd, fdp1] = fib_pair_mod(d, mod);
    const u64 fdm1 = (fdp1 + mod - fd) % mod;

    const u64 a11 = fd;
    const u64 a10 = (fdm1 + mod - 1ULL) % mod;
    const u64 a22 = (fdp1 + mod - 1ULL) % mod;

    const int v1 = std::min(vp_limited(a11, p, cap), vp_limited(a10, p, cap));

    const u64 det_left = mul_mod(a10, a22, mod);
    const u64 det_right = mul_mod(a11, a11, mod);
    const u64 det = (det_left + mod - det_right) % mod;
    const int vdet = vp_limited(det, p, cap);

    int v2 = vdet - v1;
    if (v2 < 0) {
        v2 = 0;
    }

    const int exp_total = std::min(e, v1) + std::min(e, v2);
    return pow_u64(p, exp_total);
}

struct PrimePowerData {
    u32 modulus = 0;
    u32 prime = 0;
    int exponent = 0;
    u32 period = 0;
    std::vector<u32> period_divisors;
    std::vector<u64> fixed_counts;
};

u64 lookup_fixed_count(const PrimePowerData& data, u32 d) {
    const auto it = std::lower_bound(data.period_divisors.begin(), data.period_divisors.end(), d);
    assert(it != data.period_divisors.end() && *it == d);
    const std::size_t idx = static_cast<std::size_t>(it - data.period_divisors.begin());
    return data.fixed_counts[idx];
}

u32 lcm_u32(u32 a, u32 b) {
    return static_cast<u32>((static_cast<u64>(a) / std::gcd(a, b)) * b);
}

template <typename Func>
void for_each_divisor_from_spf(u32 n, const std::vector<int>& spf, Func&& func) {
    u32 primes[16];
    int exponents[16];
    int count = 0;
    while (n > 1U) {
        const u32 p = static_cast<u32>(spf[n]);
        int e = 0;
        do {
            n /= p;
            ++e;
        } while (n > 1U && static_cast<u32>(spf[n]) == p);
        primes[count] = p;
        exponents[count] = e;
        ++count;
    }

    auto dfs = [&](auto&& self, int idx, u32 value) -> void {
        if (idx == count) {
            func(value);
            return;
        }
        u32 pe = 1U;
        for (int i = 0; i <= exponents[idx]; ++i) {
            self(self, idx + 1, value * pe);
            pe *= primes[idx];
        }
    };
    dfs(dfs, 0, 1U);
}

class Solver947 {
public:
    explicit Solver947(u32 max_m)
        : max_m_(max_m),
          max_period_(6U * max_m),
          spf_and_primes_(build_spf_and_primes(static_cast<int>(max_period_))),
          spf_(spf_and_primes_.first),
          primes_(spf_and_primes_.second),
          period_prime_(max_m_ + 1U, 0U),
          pp_index_(max_m_ + 1U, -1),
          m2_mod_(max_period_ + 1U, 0ULL) {
        precompute_period_primes();
        precompute_prime_power_data();
        precompute_m2_mod();
    }

    u64 s_mod(u32 m) const {
        int factor_indices[8];
        int factor_count = 0;

        u32 period_m = 1;
        if (m > 1) {
            u32 x = m;
            while (x > 1) {
                const u32 p = static_cast<u32>(spf_[x]);
                u32 q = 1;
                do {
                    x /= p;
                    q *= p;
                } while (x > 1 && static_cast<u32>(spf_[x]) == p);

                const int idx = pp_index_[q];
                assert(idx >= 0);
                factor_indices[factor_count++] = idx;
                period_m = lcm_u32(period_m, pp_data_[idx].period);
            }
        }

        u64 ans = 0;
        for_each_divisor_from_spf(period_m, spf_, [&](u32 d) {
            u64 fixed = 1ULL;
            for (int k = 0; k < factor_count; ++k) {
                const int idx = factor_indices[k];
                const auto& data = pp_data_[idx];
                const u32 g = std::gcd(d, data.period);
                fixed *= lookup_fixed_count(data, g);
            }

            u64 term = fixed % MOD;
            term = mul_mod(term, (static_cast<u64>(d) * d) % MOD, MOD);
            term = mul_mod(term, m2_mod_[period_m / d], MOD);
            ans = add_mod(ans, term);
        });

        return ans;
    }

    u64 S_mod(u32 M) const {
        if (M < 100'000U) {
            return S_mod_single(M);
        }
        long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
        if (cpu_count <= 1) {
            return S_mod_single(M);
        }
        u32 thread_count = static_cast<u32>(cpu_count);
        if (thread_count > 16U) {
            thread_count = 16U;
        }
        if (thread_count > M) {
            thread_count = M;
        }
        if (thread_count <= 1U) {
            return S_mod_single(M);
        }

        struct Task {
            const Solver947* solver;
            u32 begin_m;
            u32 end_m;
            u64 partial;
        };

        auto worker = [](void* arg) -> void* {
            auto* task = static_cast<Task*>(arg);
            task->partial = 0ULL;
            for (u32 m = task->begin_m; m <= task->end_m; ++m) {
                task->partial = add_mod(task->partial, task->solver->s_mod(m));
            }
            return nullptr;
        };

        std::vector<pthread_t> threads(thread_count);
        std::vector<Task> tasks(thread_count);

        const u32 base = M / thread_count;
        const u32 rem = M % thread_count;
        u32 current = 1U;
        for (u32 i = 0; i < thread_count; ++i) {
            const u32 size = base + (i < rem ? 1U : 0U);
            tasks[i] = Task{this, current, current + size - 1U, 0ULL};
            current += size;
            const int rc = ::pthread_create(&threads[i], nullptr, worker, &tasks[i]);
            assert(rc == 0);
        }

        u64 total = 0;
        for (u32 i = 0; i < thread_count; ++i) {
            const int rc = ::pthread_join(threads[i], nullptr);
            assert(rc == 0);
            total = add_mod(total, tasks[i].partial);
        }
        return total;
    }

private:
    u64 S_mod_single(u32 M) const {
        u64 total = 0;
        for (u32 m = 1; m <= M; ++m) {
            total = add_mod(total, s_mod(m));
        }
        return total;
    }

    void precompute_period_primes() {
        for (int p : primes_) {
            if (static_cast<u32>(p) > max_m_) {
                break;
            }
            period_prime_[p] = period_prime(static_cast<u32>(p), spf_);
        }
    }

    void precompute_prime_power_data() {
        pp_data_.reserve(static_cast<std::size_t>(max_m_));

        for (int pi : primes_) {
            const u32 p = static_cast<u32>(pi);
            if (p > max_m_) {
                break;
            }

            u32 q = p;
            int e = 1;
            while (q <= max_m_) {
                PrimePowerData data;
                data.modulus = q;
                data.prime = p;
                data.exponent = e;

                if (p == 2U) {
                    if (e == 1) {
                        data.period = 3U;
                    } else if (e == 2) {
                        data.period = 6U;
                    } else {
                        data.period = 3U * (1U << (e - 1));
                    }
                } else if (p == 5U) {
                    data.period = 20U * static_cast<u32>(pow_u64(5U, e - 1));
                } else {
                    data.period = period_prime_[p] * static_cast<u32>(pow_u64(p, e - 1));
                }

                data.period_divisors = divisors_from_factorization(factorize_u32(data.period, spf_), true);
                data.fixed_counts.reserve(data.period_divisors.size());
                for (u32 d : data.period_divisors) {
                    data.fixed_counts.push_back(count_fixed_prime_power(p, e, d));
                }

                const int idx = static_cast<int>(pp_data_.size());
                pp_data_.push_back(std::move(data));
                pp_index_[q] = idx;

                if (q > max_m_ / p) {
                    break;
                }
                q *= p;
                ++e;
            }
        }
    }

    void precompute_m2_mod() {
        m2_mod_[1] = 1ULL;
        for (u32 n = 2; n <= max_period_; ++n) {
            const u32 p = static_cast<u32>(spf_[n]);
            const u32 m = n / p;
            if (m % p == 0U) {
                m2_mod_[n] = m2_mod_[m];
            } else {
                const u64 p2 = (static_cast<u64>(p) * p) % MOD;
                const u64 factor = (1ULL + MOD - p2) % MOD;
                m2_mod_[n] = mul_mod(m2_mod_[m], factor, MOD);
            }
        }
    }

    u32 max_m_;
    u32 max_period_;
    std::pair<std::vector<int>, std::vector<int>> spf_and_primes_;
    const std::vector<int>& spf_;
    const std::vector<int>& primes_;

    std::vector<u32> period_prime_;
    std::vector<PrimePowerData> pp_data_;
    std::vector<int> pp_index_;
    std::vector<u64> m2_mod_;
};

void run_validations(const Solver947& solver) {
    assert(solver.s_mod(3U) == 513ULL);
    assert(solver.s_mod(10U) == 225'820ULL);
    assert(solver.S_mod(3U) == 542ULL);
    assert(solver.S_mod(10U) == 310'897ULL);
}

}  // namespace

int main() {
    Solver947 solver(1'000'000U);
    run_validations(solver);
    std::cout << solver.S_mod(1'000'000U) << '\n';
    return 0;
}
