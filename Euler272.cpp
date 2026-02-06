#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetLimit = 100'000'000'000ULL;
constexpr int kValidationRootLimit = 2'000;
constexpr int kValidationSumLimit = 5'000'000;

constexpr u64 kMinProd4Special = 53'599ULL;   // 7*13*19*31
constexpr u64 kMinProd5Special = 1'983'709ULL; // 7*13*19*31*37
constexpr u64 kPrimeBoundDenominator = 15'561ULL; // 9*7*13*19
constexpr u64 kInfinity = std::numeric_limits<u64>::max();

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

std::vector<int> build_spf(int limit) {
    std::vector<int> spf(limit + 1, 0);
    std::vector<int> primes;
    if (limit >= 1) {
        spf[1] = 1;
    }

    for (int i = 2; i <= limit; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }

        for (const int p : primes) {
            const std::int64_t composite = static_cast<std::int64_t>(i) * p;
            if (composite > limit || p > spf[i]) {
                break;
            }
            spf[static_cast<std::size_t>(composite)] = p;
        }
    }

    return spf;
}

int special_factor_count(int n, const std::vector<int>& spf) {
    int count = 0;
    while (n > 1) {
        const int p = spf[n];
        int exponent = 0;
        while (n % p == 0) {
            n /= p;
            ++exponent;
        }

        if (p == 3) {
            if (exponent >= 2) {
                ++count;
            }
        } else if (p % 3 == 1) {
            ++count;
        }
    }
    return count;
}

u64 cube_root_count_formula(int n, const std::vector<int>& spf) {
    if (n == 1) {
        return 0;
    }

    const int special = special_factor_count(n, spf);
    u64 roots = 1;
    for (int i = 0; i < special; ++i) {
        roots *= 3;
    }
    return roots;
}

u64 cube_root_count_bruteforce(u64 n) {
    u64 count = 0;
    for (u64 x = 1; x < n; ++x) {
        const u64 x2 = static_cast<u64>((static_cast<u128>(x) * x) % n);
        const u64 x3 = static_cast<u64>((static_cast<u128>(x2) * x) % n);
        if (x3 == 1) {
            ++count;
        }
    }
    return count;
}

u128 brute_sum_for_limit(int limit) {
    const std::vector<int> spf = build_spf(limit);
    u128 sum = 0;

    for (int n = 1; n <= limit; ++n) {
        if (special_factor_count(n, spf) == 5) {
            sum += static_cast<u64>(n);
        }
    }

    return sum;
}

class Euler272Solver {
public:
    explicit Euler272Solver(u64 limit) : limit_(limit) {
        build_special_primes();
        build_min_product_table();
        build_cofactor_prefix();
    }

    u128 solve(bool allow_multithreading, unsigned requested_threads = 0) const {
        unsigned threads = requested_threads;
        if (threads == 0) {
            threads = std::thread::hardware_concurrency();
            if (threads == 0) {
                threads = 1;
            }
        }

        if (!allow_multithreading) {
            threads = 1;
        }

        const u128 case_without_9 = solve_case_without_9(threads);
        const u128 case_with_9 = solve_case_with_9(threads);
        return case_without_9 + case_with_9;
    }

private:
    u64 limit_ = 0;
    u64 max_cofactor_limit_ = 0;

    std::vector<u64> special_primes_;          // primes p == 1 (mod 3)
    std::vector<std::vector<u64>> min_prod_;   // min_prod_[r][i] = product of r consecutive primes from i
    std::vector<u64> prefix_allowed_sum_;      // F1(M): sum of allowed cofactors <= M

    void build_special_primes() {
        special_primes_.clear();

        const u64 upper_u64 = limit_ / kPrimeBoundDenominator;
        if (upper_u64 < 7) {
            return;
        }

        const int upper = static_cast<int>(upper_u64);
        std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(upper) + 1, 1);
        is_prime[0] = 0;
        is_prime[1] = 0;

        for (int p = 2; static_cast<std::int64_t>(p) * p <= upper; ++p) {
            if (!is_prime[p]) {
                continue;
            }
            for (int multiple = p * p; multiple <= upper; multiple += p) {
                is_prime[multiple] = 0;
            }
        }

        for (int p = 7; p <= upper; ++p) {
            if (is_prime[p] && p % 3 == 1) {
                special_primes_.push_back(static_cast<u64>(p));
            }
        }
    }

    void build_min_product_table() {
        const int n = static_cast<int>(special_primes_.size());
        min_prod_.assign(5, std::vector<u64>(static_cast<std::size_t>(n) + 1, kInfinity));

        for (int i = 0; i <= n; ++i) {
            min_prod_[0][i] = 1;
        }

        for (int r = 1; r <= 4; ++r) {
            min_prod_[r][n] = kInfinity;
        }

        for (int i = n - 1; i >= 0; --i) {
            for (int r = 1; r <= 4; ++r) {
                const u64 tail = min_prod_[r - 1][i + 1];
                if (tail == kInfinity || special_primes_[i] > kInfinity / tail) {
                    min_prod_[r][i] = kInfinity;
                } else {
                    min_prod_[r][i] = special_primes_[i] * tail;
                }
            }
        }
    }

    void build_cofactor_prefix() {
        const u64 max_case_without_9 = limit_ / kMinProd5Special;
        const u64 max_case_with_9 = limit_ / (9ULL * kMinProd4Special);
        max_cofactor_limit_ = std::max(max_case_without_9, max_case_with_9);

        prefix_allowed_sum_.assign(static_cast<std::size_t>(max_cofactor_limit_) + 1, 0);
        if (max_cofactor_limit_ == 0) {
            return;
        }

        const int max_limit = static_cast<int>(max_cofactor_limit_);
        const std::vector<int> spf = build_spf(max_limit);

        std::vector<std::uint8_t> valid(static_cast<std::size_t>(max_limit) + 1, 0);
        valid[1] = 1;

        for (int n = 2; n <= max_limit; ++n) {
            const int p = spf[n];
            valid[n] = static_cast<std::uint8_t>(valid[n / p] && (p % 3 == 2));
        }

        for (u64 n = 1; n <= max_cofactor_limit_; ++n) {
            prefix_allowed_sum_[n] = prefix_allowed_sum_[n - 1] + (valid[n] ? n : 0);
        }
    }

    u64 sum_allowed_without_3(u64 m) const {
        if (m > max_cofactor_limit_) {
            m = max_cofactor_limit_;
        }
        return prefix_allowed_sum_[m];
    }

    u128 sum_allowed_with_optional_single_3(u64 m) const {
        if (m > max_cofactor_limit_) {
            m = max_cofactor_limit_;
        }

        const u64 no_three = sum_allowed_without_3(m);
        const u64 one_three = 3ULL * sum_allowed_without_3(m / 3);
        return static_cast<u128>(no_three) + one_three;
    }

    u128 dfs_case_without_9(int start, int remaining, u64 product) const {
        if (remaining == 0) {
            const u64 cofactor_limit = limit_ / product;
            return static_cast<u128>(product) * sum_allowed_with_optional_single_3(cofactor_limit);
        }

        const int n = static_cast<int>(special_primes_.size());
        u128 total = 0;

        for (int i = start; i < n; ++i) {
            const u64 min_rest = min_prod_[remaining - 1][i + 1];
            if (min_rest == kInfinity) {
                break;
            }
            if (product > limit_ / min_rest) {
                break;
            }

            const u64 max_prime_power = (limit_ / product) / min_rest;
            const u64 p = special_primes_[i];
            if (p > max_prime_power) {
                break;
            }

            for (u64 power = p; power <= max_prime_power;) {
                total += dfs_case_without_9(i + 1, remaining - 1, product * power);
                if (power > max_prime_power / p) {
                    break;
                }
                power *= p;
            }
        }

        return total;
    }

    u128 dfs_case_with_9(int start, int remaining, u64 product, u64 special_limit) const {
        if (remaining == 0) {
            u128 total = 0;
            const u64 max_three_power = limit_ / product;

            for (u64 three_power = 9; three_power <= max_three_power;) {
                const u64 cofactor_limit = limit_ / (product * three_power);
                total += static_cast<u128>(product) * three_power *
                         sum_allowed_without_3(cofactor_limit);
                if (three_power > max_three_power / 3) {
                    break;
                }
                three_power *= 3;
            }
            return total;
        }

        const int n = static_cast<int>(special_primes_.size());
        u128 total = 0;

        for (int i = start; i < n; ++i) {
            const u64 min_rest = min_prod_[remaining - 1][i + 1];
            if (min_rest == kInfinity) {
                break;
            }
            if (product > special_limit / min_rest) {
                break;
            }

            const u64 max_prime_power = (special_limit / product) / min_rest;
            const u64 p = special_primes_[i];
            if (p > max_prime_power) {
                break;
            }

            for (u64 power = p; power <= max_prime_power;) {
                total += dfs_case_with_9(i + 1, remaining - 1, product * power, special_limit);
                if (power > max_prime_power / p) {
                    break;
                }
                power *= p;
            }
        }

        return total;
    }

    u128 solve_case_without_9(unsigned threads) const {
        if (special_primes_.empty()) {
            return 0;
        }

        std::vector<int> roots;
        roots.reserve(special_primes_.size());

        for (int i = 0; i < static_cast<int>(special_primes_.size()); ++i) {
            const u64 min_rest = min_prod_[4][i + 1];
            if (min_rest == kInfinity) {
                break;
            }
            if (special_primes_[i] > limit_ / min_rest) {
                break;
            }
            roots.push_back(i);
        }

        if (roots.empty()) {
            return 0;
        }

        auto solve_root = [&](int root_index) {
            const u64 p = special_primes_[root_index];
            const u64 min_rest = min_prod_[4][root_index + 1];
            const u64 max_prime_power = limit_ / min_rest;

            u128 local = 0;
            for (u64 power = p; power <= max_prime_power;) {
                local += dfs_case_without_9(root_index + 1, 4, power);
                if (power > max_prime_power / p) {
                    break;
                }
                power *= p;
            }
            return local;
        };

        if (threads <= 1 || roots.size() < 2) {
            u128 total = 0;
            for (const int root : roots) {
                total += solve_root(root);
            }
            return total;
        }

        threads = std::min<unsigned>(threads, static_cast<unsigned>(roots.size()));

        std::atomic<std::size_t> next{0};
        std::vector<u128> partial(threads, 0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                u128 local = 0;
                while (true) {
                    const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= roots.size()) {
                        break;
                    }
                    local += solve_root(roots[index]);
                }
                partial[tid] = local;
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        u128 total = 0;
        for (const u128 value : partial) {
            total += value;
        }
        return total;
    }

    u128 solve_case_with_9(unsigned threads) const {
        if (limit_ < 9 || special_primes_.empty()) {
            return 0;
        }

        const u64 special_limit = limit_ / 9;

        std::vector<int> roots;
        roots.reserve(special_primes_.size());

        for (int i = 0; i < static_cast<int>(special_primes_.size()); ++i) {
            const u64 min_rest = min_prod_[3][i + 1];
            if (min_rest == kInfinity) {
                break;
            }
            if (special_primes_[i] > special_limit / min_rest) {
                break;
            }
            roots.push_back(i);
        }

        if (roots.empty()) {
            return 0;
        }

        auto solve_root = [&](int root_index) {
            const u64 p = special_primes_[root_index];
            const u64 min_rest = min_prod_[3][root_index + 1];
            const u64 max_prime_power = special_limit / min_rest;

            u128 local = 0;
            for (u64 power = p; power <= max_prime_power;) {
                local += dfs_case_with_9(root_index + 1, 3, power, special_limit);
                if (power > max_prime_power / p) {
                    break;
                }
                power *= p;
            }
            return local;
        };

        if (threads <= 1 || roots.size() < 2) {
            u128 total = 0;
            for (const int root : roots) {
                total += solve_root(root);
            }
            return total;
        }

        threads = std::min<unsigned>(threads, static_cast<unsigned>(roots.size()));

        std::atomic<std::size_t> next{0};
        std::vector<u128> partial(threads, 0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                u128 local = 0;
                while (true) {
                    const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= roots.size()) {
                        break;
                    }
                    local += solve_root(roots[index]);
                }
                partial[tid] = local;
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        u128 total = 0;
        for (const u128 value : partial) {
            total += value;
        }
        return total;
    }
};

bool run_validation_checkpoints() {
    {
        const u64 roots_mod_91 = cube_root_count_bruteforce(91);
        if (roots_mod_91 != 9) {
            std::cerr << "Checkpoint failed for n=91: got C(91)=" << (roots_mod_91 - 1)
                      << ", expected 8\n";
            return false;
        }
    }

    {
        const std::vector<int> spf = build_spf(kValidationRootLimit);
        for (int n = 1; n <= kValidationRootLimit; ++n) {
            const u64 formula = cube_root_count_formula(n, spf);
            const u64 brute = cube_root_count_bruteforce(static_cast<u64>(n));
            if (formula != brute) {
                std::cerr << "Root-count checkpoint failed for n=" << n << ": got "
                          << formula << ", expected " << brute << "\n";
                return false;
            }
        }
    }

    {
        const u128 brute = brute_sum_for_limit(kValidationSumLimit);
        const Euler272Solver solver(kValidationSumLimit);
        const u128 fast = solver.solve(false, 1);

        if (fast != brute) {
            std::cerr << "Sum checkpoint failed for n<=" << kValidationSumLimit << ": got "
                      << to_string_u128(fast) << ", expected " << to_string_u128(brute)
                      << "\n";
            return false;
        }
    }

    return true;
}

} // namespace

int main() {
    if (!run_validation_checkpoints()) {
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    const Euler272Solver solver(kTargetLimit);
    const u128 answer = solver.solve(true, threads);

    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
