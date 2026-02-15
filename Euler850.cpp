#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <vector>

using u32 = uint32_t;
using u64 = uint64_t;
using i64 = int64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

namespace {

constexpr u64 MOD = 977676779ULL;

u64 add_mod(u64 a, u64 b, u64 mod) {
    a += b;
    if (a >= mod) {
        a -= mod;
    }
    return a;
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 norm_i128_mod(i128 x, u64 mod) {
    i128 r = x % static_cast<i128>(mod);
    if (r < 0) {
        r += static_cast<i128>(mod);
    }
    return static_cast<u64>(r);
}

u64 isqrt_u64(u64 x) {
    u64 lo = 0;
    u64 hi = (x < (1ULL << 32)) ? (1ULL << 32) : x;
    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;
        if (mid <= x / mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

std::vector<u32> primes_up_to(int n) {
    std::vector<int> lp(n + 1, 0);
    std::vector<u32> primes;
    primes.reserve(static_cast<size_t>(n / 10));
    for (int i = 2; i <= n; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes.push_back(static_cast<u32>(i));
        }
        for (u32 p : primes) {
            i64 v = static_cast<i64>(i) * static_cast<i64>(p);
            if (v > n || static_cast<int>(p) > lp[i]) {
                break;
            }
            lp[static_cast<size_t>(v)] = static_cast<int>(p);
        }
    }
    return primes;
}

u64 pow_int(u64 base, int exp) {
    u64 res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= base;
    }
    return res;
}

u64 brute_T_small(int N) {
    std::vector<int> spf(N + 1, 0);
    for (int i = 2; i <= N; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            if (static_cast<i64>(i) * i <= N) {
                for (int j = i * i; j <= N; j += i) {
                    if (spf[j] == 0) {
                        spf[j] = i;
                    }
                }
            }
        }
    }

    std::vector<std::vector<std::pair<int, int>>> factors(N + 1);
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

    u64 T = 0;
    for (int k = 1; k <= N; k += 2) {
        for (int n = 1; n <= N; ++n) {
            u64 m = 1;
            for (const auto& pe : factors[n]) {
                int c = (pe.second + k - 1) / k;
                m *= pow_int(static_cast<u64>(pe.first), c);
            }
            T += static_cast<u64>(n - n / static_cast<int>(m));
        }
    }
    return T;
}

class Solver850 {
public:
    Solver850(u64 n, u64 mod, bool use_threads)
        : N_(n),
          mod2_(2 * mod),
          odd_half_((n + 1) / 2),
          use_threads_(use_threads) {
        const int lim = static_cast<int>(isqrt_u64(N_)) + 1;
        primes_ = primes_up_to(lim);
        prime_limit_ = 0;
        while (prime_limit_ < static_cast<int>(primes_.size()) &&
               static_cast<u64>(primes_[prime_limit_]) * static_cast<u64>(primes_[prime_limit_]) <= N_) {
            ++prime_limit_;
        }
    }

    u64 compute_T_mod() const {
        std::vector<u32> ps;
        std::vector<int> es;
        u64 h_mod = contribute(1ULL, 1ULL, ps, es, -1, 1ULL);
        h_mod = add_mod(h_mod, compute_branches_mod(), mod2_);

        const u64 sum_n_mod = static_cast<u64>((static_cast<u128>(N_) * (N_ + 1) / 2) % mod2_);
        const u64 a_mod = mul_mod(odd_half_ % mod2_, sum_n_mod, mod2_);
        if (a_mod >= h_mod) {
            return a_mod - h_mod;
        }
        return a_mod + mod2_ - h_mod;
    }

private:
    u64 contribute(u64 prd,
                   u64 rad,
                   const std::vector<u32>& ps,
                   const std::vector<int>& es,
                   int max_e,
                   u64 phi_mod) const {
        const u64 Np = N_ / prd;
        const int K = max_e + 2;
        const i64 shift = static_cast<i64>(K / 2);
        const i128 base = static_cast<i128>(Np / rad) * static_cast<i128>(static_cast<i64>(odd_half_) - shift);
        u64 ret = norm_i128_mod(base, mod2_);

        for (int k = K - (K & 1) - 2; k > 0; k -= 2) {
            u64 newr = 1;
            bool over = false;
            for (size_t i = 0; i < ps.size(); ++i) {
                const u64 p = ps[i];
                const int need = es[i] / k + 1;
                for (int t = 0; t < need; ++t) {
                    if (newr > Np / p) {
                        over = true;
                        break;
                    }
                    newr *= p;
                }
                if (over) {
                    break;
                }
            }
            if (over) {
                break;
            }
            ret += Np / newr;
            ret %= mod2_;
        }

        return mul_mod(ret, phi_mod, mod2_);
    }

    u64 rec(int ind,
            u64 prd,
            u64 rad,
            std::vector<u32>& ps,
            std::vector<int>& es,
            int max_e,
            u64 phi_mod) const {
        u64 ret = contribute(prd, rad, ps, es, max_e, phi_mod);

        for (int i = ind; i < prime_limit_; ++i) {
            const u64 p = primes_[i];
            if (prd > N_ / p || rad > N_ / p) {
                break;
            }
            const u64 nxt_prd = prd * p;
            const u64 nxt_rad = rad * p;
            const u64 prd_lim = N_ / nxt_rad;
            if (nxt_prd > prd_lim) {
                break;
            }

            ps.push_back(static_cast<u32>(p));
            es.push_back(-1);

            u64 cur_prd = nxt_prd;
            u64 phi_e_mod = phi_mod;
            while (cur_prd <= prd_lim) {
                const int e = es.back() + 1;
                es.back() = e;
                if (e == 0) {
                    phi_e_mod = mul_mod(phi_mod, (p - 1) % mod2_, mod2_);
                } else {
                    phi_e_mod = mul_mod(phi_e_mod, p % mod2_, mod2_);
                }

                const int child_max = std::max(max_e, e);
                ret = add_mod(ret, rec(i + 1, cur_prd, nxt_rad, ps, es, child_max, phi_e_mod), mod2_);

                if (cur_prd > prd_lim / p) {
                    break;
                }
                cur_prd *= p;
            }

            ps.pop_back();
            es.pop_back();
        }

        return ret;
    }

    u64 branch_from_prime(int idx) const {
        const u64 p = primes_[idx];
        if (p > N_ / p) {
            return 0;
        }

        const u64 prd_lim = N_ / p;
        if (p > prd_lim) {
            return 0;
        }

        std::vector<u32> ps(1, static_cast<u32>(p));
        std::vector<int> es(1, -1);
        u64 ret = 0;
        u64 cur_prd = p;
        u64 phi_e_mod = 1;
        while (cur_prd <= prd_lim) {
            const int e = es[0] + 1;
            es[0] = e;
            if (e == 0) {
                phi_e_mod = (p - 1) % mod2_;
            } else {
                phi_e_mod = mul_mod(phi_e_mod, p % mod2_, mod2_);
            }
            ret = add_mod(ret, rec(idx + 1, cur_prd, p, ps, es, e, phi_e_mod), mod2_);
            if (cur_prd > prd_lim / p) {
                break;
            }
            cur_prd *= p;
        }
        return ret;
    }

    struct WorkerArg {
        const Solver850* self = nullptr;
        std::atomic<int>* next_idx = nullptr;
        u64 partial = 0;
    };

    static void* worker_entry(void* raw) {
        auto* arg = static_cast<WorkerArg*>(raw);
        u64 local = 0;
        while (true) {
            const int idx = arg->next_idx->fetch_add(1, std::memory_order_relaxed);
            if (idx >= arg->self->prime_limit_) {
                break;
            }
            local = add_mod(local, arg->self->branch_from_prime(idx), arg->self->mod2_);
        }
        arg->partial = local;
        return nullptr;
    }

    u64 compute_branches_mod() const {
        if (prime_limit_ == 0) {
            return 0;
        }

        long cpu = ::sysconf(_SC_NPROCESSORS_ONLN);
        int thread_count = (cpu > 1) ? static_cast<int>(cpu) : 1;
        if (!use_threads_ || thread_count <= 1) {
            u64 total = 0;
            for (int i = 0; i < prime_limit_; ++i) {
                total = add_mod(total, branch_from_prime(i), mod2_);
            }
            return total;
        }

        if (thread_count > 8) {
            thread_count = 8;
        }
        if (thread_count > prime_limit_) {
            thread_count = prime_limit_;
        }

        std::atomic<int> next_idx(0);
        std::vector<pthread_t> tids(static_cast<size_t>(thread_count));
        std::vector<WorkerArg> args(static_cast<size_t>(thread_count));

        for (int t = 0; t < thread_count; ++t) {
            args[static_cast<size_t>(t)].self = this;
            args[static_cast<size_t>(t)].next_idx = &next_idx;
            const int rc = pthread_create(&tids[static_cast<size_t>(t)], nullptr, worker_entry, &args[static_cast<size_t>(t)]);
            assert(rc == 0);
        }

        u64 total = 0;
        for (int t = 0; t < thread_count; ++t) {
            const int rc = pthread_join(tids[static_cast<size_t>(t)], nullptr);
            assert(rc == 0);
            total = add_mod(total, args[static_cast<size_t>(t)].partial, mod2_);
        }
        return total;
    }

    u64 N_;
    u64 mod2_;
    u64 odd_half_;
    bool use_threads_;
    std::vector<u32> primes_;
    int prime_limit_;
};

void validate() {
    assert(brute_T_small(10) == 201ULL);
    assert(brute_T_small(30) == Solver850(30ULL, MOD, false).compute_T_mod());
    assert(Solver850(10ULL, MOD, false).compute_T_mod() == 201ULL);
    assert(Solver850(1000ULL, MOD, false).compute_T_mod() == 247375608ULL);
}

}  // namespace

int main() {
    validate();

    constexpr u64 N = 33557799775533ULL;
    const u64 t_mod = Solver850(N, MOD, true).compute_T_mod();
    const u64 answer = (t_mod - (t_mod & 1ULL)) / 2ULL;
    std::cout << (answer % MOD) << '\n';
    return 0;
}
