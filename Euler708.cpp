#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

std::vector<int> sieve_primes(const int limit) {
    if (limit < 2) {
        return {};
    }

    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;

    for (int i = 2; static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(limit); ++i) {
        if (is_prime[static_cast<std::size_t>(i)] == 0U) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0U;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(
        static_cast<long double>(limit) / std::max(1.0L, std::log(static_cast<long double>(limit)))));
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[static_cast<std::size_t>(i)] != 0U) {
            primes.push_back(i);
        }
    }
    return primes;
}

bool is_prime_trial(const u64 x, const std::vector<int>& primes) {
    if (x < 2ULL) {
        return false;
    }
    for (const int p : primes) {
        const u64 pu = static_cast<u64>(p);
        if (pu * pu > x) {
            break;
        }
        if (x % pu == 0ULL) {
            return x == pu;
        }
    }
    return true;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        out.push_back(static_cast<char>('0' + static_cast<int>(value % 10)));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

class PiSummatory {
  public:
    explicit PiSummatory(const u64 n) : n_(n), sqrt_n_(isqrt_u64(n_)), primes_(sieve_primes(static_cast<int>(sqrt_n_))) {
        next_prime_after_sqrt_ = static_cast<u64>(sqrt_n_) + 1ULL;
        while (!is_prime_trial(next_prime_after_sqrt_, primes_)) {
            ++next_prime_after_sqrt_;
        }

        for (u64 l = 1; l <= n_;) {
            const u64 w = n_ / l;
            values_.push_back(w);
            l = n_ / w + 1ULL;
        }

        const int m = static_cast<int>(values_.size());
        pi_values_.assign(static_cast<std::size_t>(m), 0ULL);
        idx_small_.assign(static_cast<std::size_t>(sqrt_n_ + 1ULL), -1);
        idx_large_.assign(static_cast<std::size_t>(sqrt_n_ + 1ULL), -1);

        for (int i = 0; i < m; ++i) {
            const u64 w = values_[static_cast<std::size_t>(i)];
            pi_values_[static_cast<std::size_t>(i)] = w - 1ULL;
            if (w <= sqrt_n_) {
                idx_small_[static_cast<std::size_t>(w)] = i;
            } else {
                idx_large_[static_cast<std::size_t>(n_ / w)] = i;
            }
        }

        for (std::size_t i = 0; i < primes_.size(); ++i) {
            const u64 p = static_cast<u64>(primes_[i]);
            const u64 p2 = p * p;
            if (p2 > n_) {
                break;
            }
            for (int j = 0; j < m && values_[static_cast<std::size_t>(j)] >= p2; ++j) {
                const u64 w = values_[static_cast<std::size_t>(j)] / p;
                const int idx = id_of(w);
                pi_values_[static_cast<std::size_t>(j)] -= pi_values_[static_cast<std::size_t>(idx)] - static_cast<u64>(i);
            }
        }
    }

    u64 n() const { return n_; }
    u64 next_prime_after_sqrt() const { return next_prime_after_sqrt_; }
    const std::vector<int>& primes() const { return primes_; }

    u64 pi(const u64 x) const {
        if (x < 2ULL) {
            return 0ULL;
        }
        return pi_values_[static_cast<std::size_t>(id_of(x))];
    }

  private:
    u64 n_ = 0ULL;
    u64 sqrt_n_ = 0ULL;
    u64 next_prime_after_sqrt_ = 0ULL;
    std::vector<u64> values_;
    std::vector<u64> pi_values_;
    std::vector<int> idx_small_;
    std::vector<int> idx_large_;
    std::vector<int> primes_;

    int id_of(const u64 x) const {
        if (x <= sqrt_n_) {
            return idx_small_[static_cast<std::size_t>(x)];
        }
        return idx_large_[static_cast<std::size_t>(n_ / x)];
    }
};

class Solver {
  public:
    explicit Solver(const u64 n) : pi_data_(n) {
        const std::size_t np = pi_data_.primes().size();
        prefix_two_.assign(np + 1ULL, 0ULL);
        for (std::size_t i = 0; i < np; ++i) {
            prefix_two_[i + 1ULL] = prefix_two_[i] + 2ULL;
        }
    }

    u128 solve(const unsigned threads) const {
        if (threads <= 1U) {
            return solve_single();
        }
        return solve_parallel(threads);
    }

  private:
    PiSummatory pi_data_;
    std::vector<u64> prefix_two_;

    u128 dfs(const u64 n, const int start) const {
        if (n < 2ULL) {
            return 0;
        }

        const int np = static_cast<int>(pi_data_.primes().size());
        if (start < np) {
            if (static_cast<u64>(pi_data_.primes()[static_cast<std::size_t>(start)]) > n) {
                return 0;
            }
        } else if (pi_data_.next_prime_after_sqrt() > n) {
            return 0;
        }

        u128 ans = static_cast<u128>(2ULL * pi_data_.pi(n) - prefix_two_[static_cast<std::size_t>(start)]);

        for (int i = start; i < np; ++i) {
            const u64 p = static_cast<u64>(pi_data_.primes()[static_cast<std::size_t>(i)]);
            if (p * p > n) {
                break;
            }
            u64 pe = p;
            u128 weight = 2;
            int exp = 1;
            while (pe <= n) {
                ans += weight * dfs(n / pe, i + 1);
                if (exp >= 2) {
                    ans += weight;
                }
                if (pe > n / p) {
                    break;
                }
                pe *= p;
                weight *= 2;
                ++exp;
            }
        }
        return ans;
    }

    u128 solve_single() const {
        return static_cast<u128>(1ULL) + dfs(pi_data_.n(), 0);
    }

    u128 solve_parallel(const unsigned threads) const {
        const auto& primes = pi_data_.primes();
        const int np = static_cast<int>(primes.size());

        int top_count = 0;
        while (top_count < np &&
               static_cast<u64>(primes[static_cast<std::size_t>(top_count)]) *
                       static_cast<u64>(primes[static_cast<std::size_t>(top_count)]) <=
                   pi_data_.n()) {
            ++top_count;
        }

        const u128 base = static_cast<u128>(2ULL * pi_data_.pi(pi_data_.n()));
        std::atomic<int> next_idx(0);
        std::vector<u128> partials(static_cast<std::size_t>(threads), 0);
        std::vector<std::thread> pool;
        pool.reserve(static_cast<std::size_t>(threads));

        for (unsigned t = 0; t < threads; ++t) {
            pool.emplace_back([&, t]() {
                u128 local = 0;
                while (true) {
                    const int i = next_idx.fetch_add(1, std::memory_order_relaxed);
                    if (i >= top_count) {
                        break;
                    }
                    const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(i)]);
                    u64 pe = p;
                    u128 weight = 2;
                    int exp = 1;
                    while (pe <= pi_data_.n()) {
                        local += weight * dfs(pi_data_.n() / pe, i + 1);
                        if (exp >= 2) {
                            local += weight;
                        }
                        if (pe > pi_data_.n() / p) {
                            break;
                        }
                        pe *= p;
                        weight *= 2;
                        ++exp;
                    }
                }
                partials[static_cast<std::size_t>(t)] = local;
            });
        }

        for (auto& th : pool) {
            th.join();
        }

        u128 ans = static_cast<u128>(1ULL) + base;
        for (const u128 v : partials) {
            ans += v;
        }
        return ans;
    }
};

u64 brute_sum(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        spf[static_cast<std::size_t>(i)] = i;
        if (static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(n)) {
            for (int j = i * i; j <= n; j += i) {
                if (spf[static_cast<std::size_t>(j)] == 0) {
                    spf[static_cast<std::size_t>(j)] = i;
                }
            }
        }
    }

    u64 sum = 1ULL;
    for (int x = 2; x <= n; ++x) {
        int t = x;
        u64 term = 1ULL;
        while (t > 1) {
            const int p = spf[static_cast<std::size_t>(t)];
            while (t % p == 0) {
                t /= p;
                term <<= 1ULL;
            }
        }
        sum += term;
    }
    return sum;
}

u128 solve_n(const u64 n, const unsigned threads) {
    const Solver solver(n);
    return solver.solve(threads);
}

void run_validations() {
    assert(static_cast<u64>(solve_n(100ULL, 1U)) == brute_sum(100));
    assert(static_cast<u64>(solve_n(50'000ULL, 1U)) == brute_sum(50'000));
    assert(static_cast<u64>(solve_n(100'000'000ULL, 1U)) == 9'613'563'919ULL);
}

}  // namespace

int main() {
    run_validations();

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0U) {
        threads = 4U;
    }

    constexpr u64 kN = 100'000'000'000'000ULL;
    const u128 ans = solve_n(kN, threads);
    std::cout << to_string_u128(ans) << '\n';
    return 0;
}
