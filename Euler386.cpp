#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using i64 = long long;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u8 = std::uint8_t;

struct Options {
    int limit = 100000000;
    int block_size = 1 << 20;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--block-size=", options.block_size)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1 && options.block_size >= 1024;
}

std::vector<int> primes_up_to(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }

    for (int p = 2; 1LL * p * p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int x = p * p; x <= limit; x += p) {
            is_prime[static_cast<std::size_t>(x)] = false;
        }
    }

    std::vector<int> primes;
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            primes.push_back(p);
        }
    }
    return primes;
}

u64 encode_pattern(const u8* exponents, const int count) {
    u64 key = static_cast<u64>(count);
    for (int i = 0; i < count; ++i) {
        key = (key << 6) | static_cast<u64>(exponents[i]);
    }
    return key;
}

u32 antichain_width_from_pattern(const u8* exponents, const int count) {
    int total_degree = 0;
    std::vector<u32> coeff(1, 1U);
    for (int idx = 0; idx < count; ++idx) {
        const int e = static_cast<int>(exponents[idx]);
        std::vector<u32> next(static_cast<std::size_t>(total_degree + e + 1), 0U);
        for (int s = 0; s <= total_degree; ++s) {
            const u32 base = coeff[static_cast<std::size_t>(s)];
            for (int add = 0; add <= e; ++add) {
                next[static_cast<std::size_t>(s + add)] += base;
            }
        }
        total_degree += e;
        coeff.swap(next);
    }
    return *std::max_element(coeff.begin(), coeff.end());
}

u32 width_from_n_by_trial_division(int n) {
    if (n == 1) {
        return 1U;
    }
    u8 exponents[16];
    int count = 0;
    int x = n;
    for (int p = 2; 1LL * p * p <= x; ++p) {
        if (x % p != 0) {
            continue;
        }
        u8 e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        exponents[count++] = e;
    }
    if (x > 1) {
        exponents[count++] = 1;
    }
    std::sort(exponents, exponents + count, std::greater<u8>());
    return antichain_width_from_pattern(exponents, count);
}

u32 brute_antichain_width(const int n) {
    std::vector<int> divisors;
    for (int d = 1; 1LL * d * d <= n; ++d) {
        if (n % d != 0) {
            continue;
        }
        divisors.push_back(d);
        if (d * d != n) {
            divisors.push_back(n / d);
        }
    }
    std::sort(divisors.begin(), divisors.end());
    const int dcnt = static_cast<int>(divisors.size());
    u32 best = 0;
    const int masks = 1 << dcnt;
    for (int mask = 1; mask < masks; ++mask) {
        const u32 bits = static_cast<u32>(__builtin_popcount(static_cast<unsigned>(mask)));
        if (bits <= best) {
            continue;
        }
        bool ok = true;
        for (int i = 0; i < dcnt && ok; ++i) {
            if (((mask >> i) & 1) == 0) {
                continue;
            }
            for (int j = i + 1; j < dcnt; ++j) {
                if (((mask >> j) & 1) == 0) {
                    continue;
                }
                const int a = divisors[static_cast<std::size_t>(i)];
                const int b = divisors[static_cast<std::size_t>(j)];
                if (a % b == 0 || b % a == 0) {
                    ok = false;
                    break;
                }
            }
        }
        if (ok) {
            best = bits;
        }
    }
    return best;
}

class Solver {
   public:
    i64 solve(const int limit, const int block_size) {
        if (limit < 1) {
            return 0;
        }

        const int root = static_cast<int>(std::sqrt(static_cast<long double>(limit)));
        const std::vector<int> primes = primes_up_to(root);

        constexpr int kMaxDistinctPrimeFactors = 10;
        std::vector<u32> remaining(static_cast<std::size_t>(block_size), 0U);
        std::vector<u8> factor_count(static_cast<std::size_t>(block_size), 0U);
        std::vector<u8> exponents(static_cast<std::size_t>(block_size * kMaxDistinctPrimeFactors), 0U);

        i64 answer = 0;
        width_cache_.clear();
        width_cache_.reserve(4096);

        for (int block_start = 1; block_start <= limit; block_start += block_size) {
            const int block_end = std::min(limit, block_start + block_size - 1);
            const int current_size = block_end - block_start + 1;

            for (int i = 0; i < current_size; ++i) {
                remaining[static_cast<std::size_t>(i)] = static_cast<u32>(block_start + i);
                factor_count[static_cast<std::size_t>(i)] = 0U;
            }

            const int local_root = static_cast<int>(std::sqrt(static_cast<long double>(block_end)));
            const auto root_it = std::upper_bound(primes.begin(), primes.end(), local_root);
            const int prime_limit = static_cast<int>(root_it - primes.begin());

            for (int pi = 0; pi < prime_limit; ++pi) {
                const int p = primes[static_cast<std::size_t>(pi)];
                i64 first = (static_cast<i64>(block_start) + p - 1) / p;
                first *= p;

                for (i64 x = first; x <= block_end; x += p) {
                    const int idx = static_cast<int>(x - block_start);
                    u32& value = remaining[static_cast<std::size_t>(idx)];
                    if (value % static_cast<u32>(p) != 0U) {
                        continue;
                    }

                    u8 exponent = 0;
                    while (value % static_cast<u32>(p) == 0U) {
                        value /= static_cast<u32>(p);
                        ++exponent;
                    }

                    const std::size_t base = static_cast<std::size_t>(idx * kMaxDistinctPrimeFactors);
                    const u8 count = factor_count[static_cast<std::size_t>(idx)];
                    exponents[base + count] = exponent;
                    factor_count[static_cast<std::size_t>(idx)] = static_cast<u8>(count + 1U);
                }
            }

            for (int idx = 0; idx < current_size; ++idx) {
                const std::size_t bidx = static_cast<std::size_t>(idx);
                const std::size_t base = static_cast<std::size_t>(idx * kMaxDistinctPrimeFactors);

                if (remaining[bidx] > 1U) {
                    const u8 count = factor_count[bidx];
                    exponents[base + count] = 1U;
                    factor_count[bidx] = static_cast<u8>(count + 1U);
                }

                const int count = static_cast<int>(factor_count[bidx]);
                if (count == 0) {
                    answer += 1;
                    continue;
                }

                std::sort(exponents.begin() + static_cast<i64>(base),
                          exponents.begin() + static_cast<i64>(base + count),
                          std::greater<u8>());

                const u64 key = encode_pattern(exponents.data() + base, count);
                const auto it = width_cache_.find(key);
                if (it != width_cache_.end()) {
                    answer += static_cast<i64>(it->second);
                    continue;
                }

                const u32 width = antichain_width_from_pattern(exponents.data() + base, count);
                width_cache_[key] = width;
                answer += static_cast<i64>(width);
            }
        }

        return answer;
    }

   private:
    std::unordered_map<u64, u32> width_cache_;
};

i64 slow_sum_by_trial_division(const int limit) {
    i64 total = 0;
    for (int n = 1; n <= limit; ++n) {
        total += static_cast<i64>(width_from_n_by_trial_division(n));
    }
    return total;
}

bool run_checkpoints() {
    if (brute_antichain_width(30) != 3U) {
        std::cerr << "Checkpoint failed: N(30)\n";
        return false;
    }

    for (int n = 1; n <= 200; ++n) {
        const u32 brute = brute_antichain_width(n);
        const u32 formula = width_from_n_by_trial_division(n);
        if (brute != formula) {
            std::cerr << "Checkpoint failed at n=" << n << " for width formula\n";
            return false;
        }
    }

    Solver solver;
    const int small_limit = 200000;
    const i64 fast = solver.solve(small_limit, 1 << 16);
    const i64 slow = slow_sum_by_trial_division(small_limit);
    if (fast != slow) {
        std::cerr << "Checkpoint failed for summatory test at limit=" << small_limit
                  << " (fast=" << fast << ", slow=" << slow << ")\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    Solver solver;
    const i64 answer = solver.solve(options.limit, options.block_size);
    std::cout << answer << '\n';
    return 0;
}
