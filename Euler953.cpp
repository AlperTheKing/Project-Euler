#include <atomic>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kInv6 = 166'666'668ULL;

u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * b) % kMod);
}

u64 sum_squares_mod(u64 n) {
    const u64 a = n % kMod;
    const u64 b = (n + 1ULL) % kMod;
    const u64 c = (2ULL * (n % kMod) + 1ULL) % kMod;
    return mod_mul(mod_mul(mod_mul(a, b), c), kInv6);
}

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 brute_sum(u64 n) {
    u64 total = 0;
    for (u64 v = 1; v <= n; ++v) {
        u64 x = v;
        u64 nim = 0;
        for (u64 p = 2; p * p <= x; ++p) {
            int parity = 0;
            while (x % p == 0) {
                x /= p;
                parity ^= 1;
            }
            if (parity != 0) {
                nim ^= p;
            }
        }
        if (x > 1) {
            nim ^= x;
        }
        if (nim == 0) {
            total += v;
        }
    }
    return total;
}

class Solver {
public:
    explicit Solver(u64 n, int threads = 0) : n_(n) {
        if (threads <= 0) {
            const unsigned hc = std::thread::hardware_concurrency();
            threads_ = static_cast<int>(hc == 0 ? 1U : hc);
        } else {
            threads_ = std::max(1, threads);
        }
        build_sieve();
    }

    u64 solve() const {
        u64 ans = 0;
        ans += solve_case(n_, 0U, false);
        if (ans >= kMod) {
            ans -= kMod;
        }
        const u64 with_two = solve_case(n_ / 2ULL, 2U, true);
        ans += with_two;
        if (ans >= kMod) {
            ans -= kMod;
        }
        return ans;
    }

private:
    u64 n_{0};
    int threads_{1};
    int sieve_limit_{0};
    std::vector<std::uint8_t> is_prime_;
    std::vector<int> odd_primes_;

    void build_sieve() {
        const u64 lim = n_ / 2ULL;
        const u64 r = isqrt_u64(lim == 0 ? 1ULL : lim);
        sieve_limit_ = static_cast<int>(2ULL * r + 100ULL);
        if (sieve_limit_ < 100) {
            sieve_limit_ = 100;
        }

        is_prime_.assign(static_cast<std::size_t>(sieve_limit_ + 1), 1U);
        is_prime_[0] = 0U;
        is_prime_[1] = 0U;

        for (int i = 2; static_cast<u64>(i) * i <= static_cast<u64>(sieve_limit_); ++i) {
            if (is_prime_[static_cast<std::size_t>(i)] == 0U) {
                continue;
            }
            for (int j = i * i; j <= sieve_limit_; j += i) {
                is_prime_[static_cast<std::size_t>(j)] = 0U;
            }
        }

        odd_primes_.clear();
        odd_primes_.reserve(static_cast<std::size_t>(sieve_limit_ / 10));
        for (int p = 3; p <= sieve_limit_; p += 2) {
            if (is_prime_[static_cast<std::size_t>(p)] != 0U) {
                odd_primes_.push_back(p);
            }
        }
    }

    bool feasible_next(u64 prod, u64 p, int rem, u64 limit) const {
        u128 v = static_cast<u128>(prod) * p;
        if (v > limit) {
            return false;
        }
        u64 t = p;
        for (int i = 0; i < rem; ++i) {
            t += 2ULL;
            v *= t;
            if (v > limit) {
                return false;
            }
        }
        return true;
    }

    u64 contribution(u64 kernel) const {
        const u64 q = n_ / kernel;
        const u64 a = isqrt_u64(q);
        const u64 ss = sum_squares_mod(a);
        return mod_mul(kernel % kMod, ss);
    }

    int max_even_odd_count(u64 limit) const {
        u64 prod = 1;
        int cnt = 0;
        for (int p : odd_primes_) {
            if (prod > limit / static_cast<u64>(p)) {
                break;
            }
            prod *= static_cast<u64>(p);
            ++cnt;
        }
        if (cnt & 1) {
            --cnt;
        }
        return std::max(0, cnt);
    }

    u64 solve_case(u64 limit, u64 target, bool include_two) const {
        if (limit == 0) {
            return 0;
        }

        u64 result = 0;
        if (!include_two && target == 0U) {
            result = contribution(1ULL);
        }

        const int max_even = max_even_odd_count(limit);
        for (int s = 2; s <= max_even; s += 2) {
            const int k = s - 1;

            if (k == 1) {
                for (int p : odd_primes_) {
                    const u64 pu = static_cast<u64>(p);
                    if (!feasible_next(1ULL, pu, k, limit)) {
                        break;
                    }
                    const u64 cand = pu ^ target;
                    if ((cand & 1ULL) == 0ULL || cand <= pu || cand > static_cast<u64>(sieve_limit_)) {
                        continue;
                    }
                    if (is_prime_[static_cast<std::size_t>(cand)] == 0U) {
                        continue;
                    }
                    if (pu > limit / cand) {
                        continue;
                    }
                    const u64 odd_kernel = pu * cand;
                    const u64 kernel = include_two ? (odd_kernel << 1ULL) : odd_kernel;
                    result += contribution(kernel);
                    if (result >= kMod) {
                        result -= kMod;
                    }
                }
                continue;
            }

            std::vector<int> roots;
            roots.reserve(1024);
            for (int i = 0; i < static_cast<int>(odd_primes_.size()); ++i) {
                const u64 p = static_cast<u64>(odd_primes_[static_cast<std::size_t>(i)]);
                if (!feasible_next(1ULL, p, k, limit)) {
                    break;
                }
                roots.push_back(i);
            }
            if (roots.empty()) {
                continue;
            }

            const int workers = std::min(threads_, static_cast<int>(roots.size()));
            std::vector<std::thread> pool;
            pool.reserve(static_cast<std::size_t>(workers));
            std::vector<u64> partial(static_cast<std::size_t>(workers), 0ULL);
            std::atomic<int> next_idx{0};

            auto dfs = [&](auto&& self,
                           int start,
                           int rem,
                           u64 prod,
                           u64 xr,
                           int last,
                           u64& local_sum) -> void {
                if (rem == 0) {
                    const u64 cand = xr ^ target;
                    if ((cand & 1ULL) == 0ULL || cand <= static_cast<u64>(last) || cand > static_cast<u64>(sieve_limit_)) {
                        return;
                    }
                    if (is_prime_[static_cast<std::size_t>(cand)] == 0U || prod > limit / cand) {
                        return;
                    }
                    const u64 odd_kernel = prod * cand;
                    const u64 kernel = include_two ? (odd_kernel << 1ULL) : odd_kernel;
                    local_sum += contribution(kernel);
                    if (local_sum >= kMod) {
                        local_sum %= kMod;
                    }
                    return;
                }

                for (int i = start; i < static_cast<int>(odd_primes_.size()); ++i) {
                    const u64 p = static_cast<u64>(odd_primes_[static_cast<std::size_t>(i)]);
                    if (!feasible_next(prod, p, rem, limit)) {
                        break;
                    }
                    self(self,
                         i + 1,
                         rem - 1,
                         prod * p,
                         xr ^ p,
                         static_cast<int>(p),
                         local_sum);
                }
            };

            for (int t = 0; t < workers; ++t) {
                pool.emplace_back([&, t]() {
                    u64 local_sum = 0;
                    while (true) {
                        const int pos = next_idx.fetch_add(1, std::memory_order_relaxed);
                        if (pos >= static_cast<int>(roots.size())) {
                            break;
                        }
                        const int ridx = roots[static_cast<std::size_t>(pos)];
                        const u64 p = static_cast<u64>(odd_primes_[static_cast<std::size_t>(ridx)]);
                        dfs(dfs,
                            ridx + 1,
                            k - 1,
                            p,
                            p,
                            static_cast<int>(p),
                            local_sum);
                    }
                    partial[static_cast<std::size_t>(t)] = local_sum % kMod;
                });
            }

            for (std::thread& th : pool) {
                th.join();
            }

            for (u64 v : partial) {
                result += v;
                result %= kMod;
            }
        }

        return result % kMod;
    }
};

void run_validations() {
    {
        const Solver solver10(10, 1);
        assert(solver10.solve() == 14ULL);
    }
    {
        const Solver solver100(100, 1);
        assert(solver100.solve() == 455ULL);
    }
    {
        constexpr u64 kCheckN = 2'000ULL;
        const Solver solver(kCheckN, 1);
        const u64 fast = solver.solve();
        const u64 slow = brute_sum(kCheckN) % kMod;
        assert(fast == slow);
    }
}

}  // namespace

int main() {
    run_validations();
    constexpr u64 kN = 100'000'000'000'000ULL;
    Solver solver(kN);
    std::cout << solver.solve() << '\n';
    return 0;
}
