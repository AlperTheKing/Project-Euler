#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct PrimeNeed {
    u64 prime;
    int need;
};

u64 pow_capped(u64 base, int exp, u64 cap) {
    u64 result = 1;
    for (int i = 0; i < exp; ++i) {
        if (result > cap / base) {
            return cap + 1;
        }
        result *= base;
    }
    return result;
}

u64 mul_pow_capped(u64 value, u64 base, int exp, u64 cap) {
    for (int i = 0; i < exp; ++i) {
        if (value > cap / base) {
            return cap + 1;
        }
        value *= base;
    }
    return value;
}

std::vector<int> sieve_primes(int max_n) {
    std::vector<bool> is_prime(max_n + 1, true);
    is_prime[0] = false;
    is_prime[1] = false;
    for (int i = 2; 1LL * i * i <= max_n; ++i) {
        if (!is_prime[i]) {
            continue;
        }
        for (int j = i * i; j <= max_n; j += i) {
            is_prime[j] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= max_n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        out.push_back(static_cast<char>('0' + (value % 10)));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

class TauSolver {
  public:
    explicit TauSolver(u64 limit) : limit_(limit), primes_(sieve_primes(300000)) {}

    std::unordered_map<u64, u64> compute_minimals() {
        best_.clear();
        factor_cache_.clear();
        std::vector<int> exponents;
        dfs_vectors(0, 64, 1, 1, exponents);
        return best_;
    }

  private:
    u64 limit_;
    std::vector<int> primes_;
    std::unordered_map<u64, std::vector<PrimeNeed>> factor_cache_;
    std::unordered_map<u64, u64> best_;

    const std::vector<PrimeNeed>& factorize(u64 x) {
        auto it = factor_cache_.find(x);
        if (it != factor_cache_.end()) {
            return it->second;
        }

        std::vector<PrimeNeed> factors;
        u64 value = x;
        for (int prime : primes_) {
            const u64 p = static_cast<u64>(prime);
            if (p * p > value) {
                break;
            }
            if (value % p != 0) {
                continue;
            }
            int exponent = 0;
            while (value % p == 0) {
                value /= p;
                ++exponent;
            }
            factors.push_back({p, exponent});
        }
        if (value > 1) {
            factors.push_back({value, 1});
        }

        auto [pos, _] = factor_cache_.emplace(x, std::move(factors));
        return pos->second;
    }

    std::vector<u64> build_extra_primes(int count, const std::vector<PrimeNeed>& mandatory) const {
        std::unordered_set<u64> blocked;
        blocked.reserve(mandatory.size() * 2 + 1);
        for (const PrimeNeed& req : mandatory) {
            blocked.insert(req.prime);
        }

        std::vector<u64> extras;
        extras.reserve(count);
        for (int prime : primes_) {
            const u64 p = static_cast<u64>(prime);
            if (blocked.find(p) != blocked.end()) {
                continue;
            }
            extras.push_back(p);
            if (static_cast<int>(extras.size()) == count) {
                break;
            }
        }
        return extras;
    }

    u64 minimal_value_for_vector(const std::vector<int>& exponents, const std::vector<PrimeNeed>& factors) {
        const int r = static_cast<int>(exponents.size());
        const int t = static_cast<int>(factors.size());
        if (t > r) {
            return limit_ + 1;
        }

        std::vector<PrimeNeed> mandatory = factors;
        std::sort(mandatory.begin(), mandatory.end(), [](const PrimeNeed& lhs, const PrimeNeed& rhs) {
            if (lhs.prime != rhs.prime) {
                return lhs.prime > rhs.prime;
            }
            return lhs.need > rhs.need;
        });

        for (const PrimeNeed& req : mandatory) {
            if (exponents.empty() || exponents.front() < req.need) {
                return limit_ + 1;
            }
        }

        const std::vector<u64> extras = build_extra_primes(r - t, mandatory);
        std::vector<bool> used(r, false);
        u64 best = limit_ + 1;

        const auto dfs_assign = [&](auto&& self, int idx, u64 current) -> void {
            if (current >= best) {
                return;
            }

            if (idx == t) {
                u64 value = current;
                int extra_idx = 0;
                for (int i = 0; i < r; ++i) {
                    if (used[i]) {
                        continue;
                    }
                    value = mul_pow_capped(value, extras[extra_idx], exponents[i], limit_);
                    if (value > limit_ || value >= best) {
                        return;
                    }
                    ++extra_idx;
                }
                best = std::min(best, value);
                return;
            }

            int previous_exp = -1;
            for (int i = 0; i < r; ++i) {
                if (used[i]) {
                    continue;
                }
                const int exp = exponents[i];
                if (exp < mandatory[idx].need) {
                    continue;
                }
                if (exp == previous_exp) {
                    continue;
                }
                previous_exp = exp;

                const u64 multiplied = mul_pow_capped(current, mandatory[idx].prime, exp, limit_);
                if (multiplied > limit_ || multiplied >= best) {
                    continue;
                }
                used[i] = true;
                self(self, idx + 1, multiplied);
                used[i] = false;
            }
        };

        dfs_assign(dfs_assign, 0, 1);
        return best;
    }

    void dfs_vectors(int prime_idx, int max_exp, u64 current_value, u64 tau_value, std::vector<int>& exponents) {
        const std::vector<PrimeNeed>& factors = factorize(tau_value);
        const u64 candidate = minimal_value_for_vector(exponents, factors);
        if (candidate <= limit_) {
            auto it = best_.find(tau_value);
            if (it == best_.end() || candidate < it->second) {
                best_[tau_value] = candidate;
            }
        }

        if (prime_idx >= static_cast<int>(primes_.size())) {
            return;
        }

        const u64 prime = static_cast<u64>(primes_[prime_idx]);
        u64 value = current_value;
        for (int exp = 1; exp <= max_exp; ++exp) {
            if (value > limit_ / prime) {
                break;
            }
            value *= prime;

            if (tau_value > std::numeric_limits<u64>::max() / static_cast<u64>(exp + 1)) {
                break;
            }
            const u64 new_tau = tau_value * static_cast<u64>(exp + 1);

            exponents.push_back(exp);
            dfs_vectors(prime_idx + 1, exp, value, new_tau, exponents);
            exponents.pop_back();
        }
    }
};

u128 sum_values(const std::unordered_map<u64, u64>& values) {
    u128 total = 0;
    for (const auto& [_, v] : values) {
        total += static_cast<u128>(v);
    }
    return total;
}

std::unordered_map<u64, u64> brute_minimals(int limit) {
    std::vector<int> tau(limit + 1, 0);
    for (int d = 1; d <= limit; ++d) {
        for (int multiple = d; multiple <= limit; multiple += d) {
            ++tau[multiple];
        }
    }

    std::unordered_map<u64, u64> best;
    for (int x = 1; x <= limit; ++x) {
        const int k = tau[x];
        if (x % k != 0) {
            continue;
        }
        const auto it = best.find(static_cast<u64>(k));
        if (it == best.end() || static_cast<u64>(x) < it->second) {
            best[static_cast<u64>(k)] = static_cast<u64>(x);
        }
    }
    return best;
}

void validate() {
    {
        TauSolver solver(1000);
        const auto minima = solver.compute_minimals();
        assert(minima.at(8) == 24);
        assert(minima.at(12) == 60);
        assert(minima.at(16) == 384);
        assert(sum_values(minima) == 3189);
    }

    {
        constexpr int kSmallLimit = 20000;
        TauSolver solver(kSmallLimit);
        const auto fast = solver.compute_minimals();
        const auto brute = brute_minimals(kSmallLimit);
        assert(fast.size() == brute.size());
        for (const auto& [k, value] : brute) {
            const auto it = fast.find(k);
            assert(it != fast.end());
            assert(it->second == value);
        }
    }
}

}  // namespace

int main() {
    validate();

    constexpr u64 kLimit = 10'000'000'000'000'000ULL;
    TauSolver solver(kLimit);
    const auto minima = solver.compute_minimals();
    std::cout << to_string_u128(sum_values(minima)) << '\n';
    return 0;
}
