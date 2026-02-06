#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultLimit = 100'000'000'000ULL;
constexpr u64 kSampleN = 10'000ULL;
constexpr u64 kSampleExpectedF = 36ULL;
constexpr u64 kGeometryCrossCheckMaxN = 400ULL;
constexpr u64 kBruteSumCrossCheckLimit = 2'000'000ULL;
constexpr u64 kThreadConsistencyLimit = 10'000'000'000ULL;
constexpr u64 kExponentOneMinOtherProduct = 21'125ULL;  // 5^3 * 13^2

struct Options {
    u64 limit = kDefaultLimit;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        u64 limit = 0;
        if (parse_u64_after_prefix(arg, "--limit=", limit)) {
            options.limit = limit;
            continue;
        }

        unsigned threads = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", threads)) {
            options.requested_threads = threads;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.limit == 0ULL) {
        std::cerr << "--limit must be >= 1.\n";
        return false;
    }

    return true;
}

u64 isqrt_u64(u64 x) {
    if (x == 0ULL) {
        return 0ULL;
    }
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 pow_with_limit(u64 base, int exponent, u64 limit) {
    u64 result = 1ULL;
    for (int i = 0; i < exponent; ++i) {
        if (result > limit / base) {
            return limit + 1ULL;
        }
        result *= base;
    }
    return result;
}

std::vector<int> sieve_primes(int limit) {
    if (limit < 2) {
        return {};
    }

    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit) + 1ULL, 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(limit)));
    for (int p = 2; p <= root; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (std::int64_t m = static_cast<std::int64_t>(p) * p; m <= limit; m += p) {
            is_prime[static_cast<std::size_t>(m)] = 0U;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(
        static_cast<long double>(limit) /
        std::max(1.0L, std::log(static_cast<long double>(limit)))));
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] != 0U) {
            primes.push_back(p);
        }
    }
    return primes;
}

std::vector<int> primes_1mod4_up_to(int limit) {
    const std::vector<int> primes = sieve_primes(limit);
    std::vector<int> out;
    out.reserve(primes.size() / 2ULL + 1ULL);
    for (const int p : primes) {
        if ((p % 4) == 1) {
            out.push_back(p);
        }
    }
    return out;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2'000ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

u64 f_value_from_factorization(u64 n, const std::vector<int>& primes) {
    u64 m = n;
    u64 product = 1ULL;

    for (const int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        if (p > m / p) {
            break;
        }
        if ((m % p) != 0ULL) {
            continue;
        }

        int exponent = 0;
        while ((m % p) == 0ULL) {
            m /= p;
            ++exponent;
        }
        if ((p & 3ULL) == 1ULL) {
            product *= static_cast<u64>(2 * exponent + 1);
        }
    }

    if (m > 1ULL && (m & 3ULL) == 1ULL) {
        product *= 3ULL;
    }

    return 4ULL * product;
}

u64 brute_circle_count(u64 n) {
    const u64 target = n * n;
    const std::int64_t nn = static_cast<std::int64_t>(n);

    u64 count = 0ULL;
    for (std::int64_t x = -nn; x <= nn; ++x) {
        const std::int64_t rem = static_cast<std::int64_t>(target) - x * x;
        if (rem < 0) {
            continue;
        }

        const u64 y = isqrt_u64(static_cast<u64>(rem));
        if (y * y != static_cast<u64>(rem)) {
            continue;
        }
        count += (y == 0ULL ? 1ULL : 2ULL);
    }

    return count;
}

std::vector<u32> build_spf(u64 max_n) {
    std::vector<u32> spf(static_cast<std::size_t>(max_n) + 1ULL, 0U);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(
        static_cast<long double>(max_n) /
        std::max(1.0L, std::log(std::max(2.0L, static_cast<long double>(max_n))))));

    for (u64 i = 2ULL; i <= max_n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0U) {
            spf[static_cast<std::size_t>(i)] = static_cast<u32>(i);
            primes.push_back(static_cast<u32>(i));
        }
        for (const u32 p : primes) {
            const u64 m = static_cast<u64>(p) * i;
            if (m > max_n) {
                break;
            }
            spf[static_cast<std::size_t>(m)] = p;
            if (p == spf[static_cast<std::size_t>(i)]) {
                break;
            }
        }
    }

    return spf;
}

u128 brute_sum_direct(u64 limit) {
    if (limit == 0ULL) {
        return 0;
    }

    const std::vector<u32> spf = build_spf(limit);
    u128 sum = 0;

    for (u64 n = 1ULL; n <= limit; ++n) {
        u64 m = n;
        u64 product = 1ULL;
        while (m > 1ULL) {
            const u32 p = spf[static_cast<std::size_t>(m)];
            int exponent = 0;
            while ((m % p) == 0ULL) {
                m /= p;
                ++exponent;
            }
            if ((p & 3U) == 1U) {
                product *= static_cast<u64>(2 * exponent + 1);
            }
        }

        const u64 f = 4ULL * product;
        if (f == 420ULL) {
            sum += n;
        }
    }

    return sum;
}

bool tail_fits(const std::vector<int>& ordered_exponents,
               std::size_t next_pos,
               std::size_t next_prime_idx,
               u64 remaining_limit,
               const std::vector<int>& primes_1mod4) {
    std::size_t idx = next_prime_idx;
    for (std::size_t pos = next_pos; pos < ordered_exponents.size(); ++pos, ++idx) {
        if (idx >= primes_1mod4.size()) {
            return false;
        }
        const u64 p = static_cast<u64>(primes_1mod4[idx]);
        const u64 pe = pow_with_limit(p, ordered_exponents[pos], remaining_limit);
        if (pe > remaining_limit) {
            return false;
        }
        remaining_limit /= pe;
    }
    return true;
}

void enumerate_core_for_order(const std::vector<int>& ordered_exponents,
                              std::size_t pos,
                              std::size_t start_prime_idx,
                              u64 current_product,
                              u64 limit,
                              const std::vector<int>& primes_1mod4,
                              std::vector<u64>& out) {
    if (pos == ordered_exponents.size()) {
        out.push_back(current_product);
        return;
    }

    const std::size_t remaining_slots = ordered_exponents.size() - pos;
    if (start_prime_idx + remaining_slots > primes_1mod4.size()) {
        return;
    }

    const u64 head_limit = limit / current_product;
    for (std::size_t i = start_prime_idx; i + remaining_slots <= primes_1mod4.size(); ++i) {
        const u64 p = static_cast<u64>(primes_1mod4[i]);
        const u64 pe = pow_with_limit(p, ordered_exponents[pos], head_limit);
        if (pe > head_limit) {
            break;
        }

        const u64 next_product = current_product * pe;
        if (pos + 1ULL < ordered_exponents.size()) {
            if (!tail_fits(ordered_exponents,
                           pos + 1ULL,
                           i + 1ULL,
                           limit / next_product,
                           primes_1mod4)) {
                break;
            }
        }

        enumerate_core_for_order(ordered_exponents,
                                 pos + 1ULL,
                                 i + 1ULL,
                                 next_product,
                                 limit,
                                 primes_1mod4,
                                 out);
    }
}

std::vector<u64> enumerate_core_numbers(u64 limit, const std::vector<int>& primes_1mod4) {
    // 4 * Π(2a_i + 1) = 420  =>  Π(2a_i + 1) = 105.
    const std::vector<std::vector<int>> exponent_patterns = {
        {52},
        {17, 1},
        {10, 2},
        {7, 3},
        {3, 2, 1},
    };

    std::vector<u64> core_numbers;
    core_numbers.reserve(200'000);

    for (const auto& pattern : exponent_patterns) {
        std::vector<int> permutation = pattern;
        std::sort(permutation.begin(), permutation.end());
        do {
            enumerate_core_for_order(permutation,
                                     0,
                                     0,
                                     1ULL,
                                     limit,
                                     primes_1mod4,
                                     core_numbers);
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }

    std::sort(core_numbers.begin(), core_numbers.end());
    core_numbers.erase(std::unique(core_numbers.begin(), core_numbers.end()), core_numbers.end());
    return core_numbers;
}

std::vector<u128> build_allowed_multiplier_prefix(u64 max_multiplier) {
    std::vector<u128> prefix(static_cast<std::size_t>(max_multiplier) + 1ULL, 0);
    if (max_multiplier == 0ULL) {
        return prefix;
    }

    const std::vector<u32> spf = build_spf(max_multiplier);
    std::vector<std::uint8_t> valid(static_cast<std::size_t>(max_multiplier) + 1ULL, 0U);
    valid[1] = 1U;

    for (u64 n = 2ULL; n <= max_multiplier; ++n) {
        const u32 p = spf[static_cast<std::size_t>(n)];
        const u64 m = n / static_cast<u64>(p);
        const bool allowed_prime = (p == 2U) || ((p & 3U) == 3U);
        valid[static_cast<std::size_t>(n)] =
            static_cast<std::uint8_t>(allowed_prime && valid[static_cast<std::size_t>(m)] != 0U);
    }

    u128 running = 0;
    for (u64 n = 1ULL; n <= max_multiplier; ++n) {
        if (valid[static_cast<std::size_t>(n)] != 0U) {
            running += n;
        }
        prefix[static_cast<std::size_t>(n)] = running;
    }
    return prefix;
}

u128 accumulate_range(const std::vector<u64>& core_numbers,
                     std::size_t begin,
                     std::size_t end,
                     u64 limit,
                     const std::vector<u128>& allowed_prefix) {
    u128 local = 0;
    for (std::size_t i = begin; i < end; ++i) {
        const u64 core = core_numbers[i];
        const u64 max_multiplier = limit / core;
        local += static_cast<u128>(core) * allowed_prefix[static_cast<std::size_t>(max_multiplier)];
    }
    return local;
}

u64 estimate_prime_bound(u64 limit) {
    const u64 by_exp_one = limit / kExponentOneMinOtherProduct + 100ULL;
    return std::max<u64>(200ULL, by_exp_one);
}

u128 solve_sum(u64 limit, bool allow_multithreading, unsigned requested_threads) {
    const u64 prime_bound_u64 = estimate_prime_bound(limit);
    if (prime_bound_u64 > static_cast<u64>(std::numeric_limits<int>::max())) {
        std::cerr << "Prime bound is too large for this build.\n";
        return 0;
    }

    const std::vector<int> primes_1mod4 = primes_1mod4_up_to(static_cast<int>(prime_bound_u64));
    const std::vector<u64> core_numbers = enumerate_core_numbers(limit, primes_1mod4);
    if (core_numbers.empty()) {
        return 0;
    }

    const u64 min_core = core_numbers.front();
    const u64 max_multiplier = limit / min_core;
    const std::vector<u128> allowed_prefix = build_allowed_multiplier_prefix(max_multiplier);

    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, core_numbers.size());
    if (threads == 1U) {
        return accumulate_range(core_numbers, 0, core_numbers.size(), limit, allowed_prefix);
    }

    std::vector<std::thread> pool;
    std::vector<u128> partial(threads, 0);
    pool.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        const std::size_t begin = core_numbers.size() * t / threads;
        const std::size_t end = core_numbers.size() * (t + 1ULL) / threads;
        pool.emplace_back([&, t, begin, end]() {
            partial[t] = accumulate_range(core_numbers, begin, end, limit, allowed_prefix);
        });
    }
    for (std::thread& th : pool) {
        th.join();
    }

    u128 total = 0;
    for (const u128 x : partial) {
        total += x;
    }
    return total;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool run_checkpoints(const Options& options) {
    const std::vector<int> factor_primes = sieve_primes(1'000'000);

    const u64 sample_f = f_value_from_factorization(kSampleN, factor_primes);
    if (sample_f != kSampleExpectedF) {
        std::cerr << "Checkpoint failed: f(" << kSampleN << ") expected " << kSampleExpectedF
                  << ", got " << sample_f << '\n';
        return false;
    }
    std::cout << "Checkpoint OK: f(" << kSampleN << ") = " << sample_f << '\n';

    for (u64 n = 1ULL; n <= kGeometryCrossCheckMaxN; ++n) {
        const u64 brute = brute_circle_count(n);
        const u64 formula = f_value_from_factorization(n, factor_primes);
        if (brute != formula) {
            std::cerr << "Checkpoint failed: geometry/formula mismatch at N=" << n
                      << ", brute=" << brute << ", formula=" << formula << '\n';
            return false;
        }
    }
    std::cout << "Checkpoint OK: geometry cross-check for N <= " << kGeometryCrossCheckMaxN
              << '\n';

    const u64 brute_limit = std::min<u64>(kBruteSumCrossCheckLimit, options.limit);
    const u128 brute_expected = brute_sum_direct(brute_limit);
    const u128 brute_fast = solve_sum(brute_limit, false, 1U);
    if (brute_expected != brute_fast) {
        std::cerr << "Checkpoint failed: brute sum mismatch at limit=" << brute_limit
                  << ", brute=" << to_string_u128(brute_expected)
                  << ", fast=" << to_string_u128(brute_fast) << '\n';
        return false;
    }
    std::cout << "Checkpoint OK: brute sum cross-check <= " << brute_limit
              << " gives " << to_string_u128(brute_fast) << '\n';

    if (options.allow_multithreading) {
        const u64 tc_limit = std::min<u64>(kThreadConsistencyLimit, options.limit);
        const u128 single = solve_sum(tc_limit, false, 1U);
        const u128 multi = solve_sum(tc_limit, true, options.requested_threads);
        if (single != multi) {
            std::cerr << "Checkpoint failed: thread consistency mismatch at limit=" << tc_limit
                      << ", single=" << to_string_u128(single)
                      << ", multi=" << to_string_u128(multi) << '\n';
            return false;
        }
        std::cout << "Checkpoint OK: threaded consistency at limit=" << tc_limit
                  << " gives " << to_string_u128(single) << '\n';
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const auto start = std::chrono::steady_clock::now();

    if (options.run_checkpoints) {
        if (!run_checkpoints(options)) {
            return 1;
        }
    }

    const u128 answer =
        solve_sum(options.limit, options.allow_multithreading, options.requested_threads);

    const auto finish = std::chrono::steady_clock::now();
    const std::chrono::duration<long double> elapsed = finish - start;

    std::cout << "Answer: " << to_string_u128(answer) << '\n';
    std::cout << "S(" << options.limit << ") = " << to_string_u128(answer) << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " s\n";

    return 0;
}
