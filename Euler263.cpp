#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultLimit = 2'000'000'000ULL;
constexpr unsigned kDefaultTargetCount = 4U;
constexpr u64 kDefaultSegmentSpan = 8'000'000ULL;

constexpr int kPracticalCrossCheckMax = 400;
constexpr u64 kKnownParadise1 = 219'869'980ULL;
constexpr u64 kKnownParadise2 = 312'501'820ULL;

struct Options {
    u64 limit = kDefaultLimit;
    unsigned target_count = kDefaultTargetCount;
    u64 segment_span = kDefaultSegmentSpan;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
};

struct SearchResult {
    std::vector<u64> paradises;
    u64 sum = 0ULL;
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

    u64 parsed = 0ULL;
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

bool parse_unsigned_after_prefix(const std::string& arg, const char* prefix, unsigned& value) {
    u64 parsed = 0ULL;
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

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--limit=", parsed_u64)) {
            options.limit = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--segment=", parsed_u64)) {
            options.segment_span = parsed_u64;
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--target=", parsed_unsigned)) {
            options.target_count = parsed_unsigned;
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.limit < 30ULL) {
        std::cerr << "--limit must be at least 30.\n";
        return false;
    }
    if (options.segment_span < 10ULL) {
        std::cerr << "--segment must be at least 10.\n";
        return false;
    }
    if (options.target_count == 0U) {
        std::cerr << "--target must be at least 1.\n";
        return false;
    }

    return true;
}

u64 isqrt_u64(u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
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
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] != 0U) {
            primes.push_back(p);
        }
    }
    return primes;
}

std::vector<std::pair<u64, int>> factorize(u64 n, const std::vector<int>& primes) {
    std::vector<std::pair<u64, int>> factors;
    if (n <= 1ULL) {
        return factors;
    }

    u64 m = n;
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
        factors.emplace_back(p, exponent);
    }

    if (m > 1ULL) {
        factors.emplace_back(m, 1);
    }

    return factors;
}

u128 sigma_prime_power(u64 prime, int exponent) {
    u128 sum = 1;
    u128 power = 1;
    for (int i = 0; i < exponent; ++i) {
        power *= static_cast<u128>(prime);
        sum += power;
    }
    return sum;
}

bool is_practical(u64 n, const std::vector<int>& factor_primes) {
    if (n == 1ULL) {
        return true;
    }
    if ((n & 1ULL) == 1ULL) {
        return false;
    }

    const std::vector<std::pair<u64, int>> factors = factorize(n, factor_primes);
    if (factors.empty() || factors[0].first != 2ULL) {
        return false;
    }

    u128 sigma_prefix = sigma_prime_power(factors[0].first, factors[0].second);
    for (std::size_t i = 1; i < factors.size(); ++i) {
        const u64 p = factors[i].first;
        if (static_cast<u128>(p) > sigma_prefix + 1U) {
            return false;
        }
        sigma_prefix *= sigma_prime_power(factors[i].first, factors[i].second);
    }

    return true;
}

void generate_divisors_rec(const std::vector<std::pair<u64, int>>& factors,
                           std::size_t index,
                           u64 current,
                           std::vector<u64>& divisors) {
    if (index == factors.size()) {
        divisors.push_back(current);
        return;
    }

    const u64 p = factors[index].first;
    const int exponent = factors[index].second;

    u64 value = 1ULL;
    for (int e = 0; e <= exponent; ++e) {
        generate_divisors_rec(factors, index + 1ULL, current * value, divisors);
        if (e != exponent) {
            value *= p;
        }
    }
}

bool is_practical_bruteforce(u64 n, const std::vector<int>& factor_primes) {
    if (n == 0ULL) {
        return false;
    }
    if (n == 1ULL) {
        return true;
    }

    const std::vector<std::pair<u64, int>> factors = factorize(n, factor_primes);
    std::vector<u64> divisors;
    divisors.reserve(128);
    generate_divisors_rec(factors, 0ULL, 1ULL, divisors);
    std::sort(divisors.begin(), divisors.end());

    std::vector<std::uint8_t> reachable(static_cast<std::size_t>(n) + 1ULL, 0U);
    reachable[0] = 1U;
    for (const u64 d : divisors) {
        for (u64 s = n; s >= d; --s) {
            if (reachable[static_cast<std::size_t>(s - d)] != 0U) {
                reachable[static_cast<std::size_t>(s)] = 1U;
            }
            if (s == d) {
                break;
            }
        }
    }

    for (u64 k = 1ULL; k <= n; ++k) {
        if (reachable[static_cast<std::size_t>(k)] == 0U) {
            return false;
        }
    }
    return true;
}

bool check_practical_criterion(const std::vector<int>& factor_primes) {
    for (int n = 1; n <= kPracticalCrossCheckMax; ++n) {
        const bool criterion = is_practical(static_cast<u64>(n), factor_primes);
        const bool brute = is_practical_bruteforce(static_cast<u64>(n), factor_primes);
        if (criterion != brute) {
            std::cerr << "Practical criterion mismatch at n=" << n
                      << " criterion=" << criterion
                      << " brute=" << brute << '\n';
            return false;
        }
    }
    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2ULL) {
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

std::vector<u64> sieve_segment(u64 low, u64 high, const std::vector<int>& base_primes) {
    if (high < 2ULL || low > high) {
        return {};
    }

    const u64 odd_low = std::max<u64>(3ULL, low | 1ULL);
    if (odd_low > high) {
        return {};
    }

    const std::size_t odd_count =
        static_cast<std::size_t>((high - odd_low) / 2ULL + 1ULL);
    std::vector<std::uint8_t> is_prime(odd_count, 1U);

    for (const int p_int : base_primes) {
        const u64 p = static_cast<u64>(p_int);
        if (p == 2ULL) {
            continue;
        }
        if (p > high / p) {
            break;
        }

        u64 start = odd_low;
        const u64 rem = start % p;
        if (rem != 0ULL) {
            start += (p - rem);
        }
        if (start < p * p) {
            start = p * p;
        }
        if ((start & 1ULL) == 0ULL) {
            start += p;
        }

        for (u64 m = start; m <= high; m += 2ULL * p) {
            const std::size_t idx = static_cast<std::size_t>((m - odd_low) / 2ULL);
            is_prime[idx] = 0U;
        }
    }

    std::vector<u64> primes;
    primes.reserve(odd_count / 8ULL + 16ULL);
    for (std::size_t i = 0; i < odd_count; ++i) {
        if (is_prime[i] != 0U) {
            primes.push_back(odd_low + 2ULL * static_cast<u64>(i));
        }
    }
    return primes;
}

bool is_engineers_paradise_candidate(u64 n, const std::vector<int>& factor_primes) {
    if (n < 9ULL) {
        return false;
    }

    static constexpr int kOffsets[5] = {-8, -4, 0, 4, 8};
    for (const int offset : kOffsets) {
        const u64 value = static_cast<u64>(static_cast<std::int64_t>(n) + offset);
        if (!is_practical(value, factor_primes)) {
            return false;
        }
    }
    return true;
}

SearchResult find_paradises(u64 limit,
                            unsigned target_count,
                            u64 segment_span,
                            bool allow_multithreading,
                            unsigned requested_threads,
                            const std::vector<int>& base_primes,
                            const std::vector<int>& factor_primes) {
    SearchResult result;
    result.paradises.reserve(target_count);

    const u64 prime_limit = limit + 9ULL;
    u64 segment_low = 3ULL;
    const std::size_t total_segments =
        static_cast<std::size_t>((prime_limit >= segment_low)
                                     ? ((prime_limit - segment_low) / segment_span + 1ULL)
                                     : 0ULL);

    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, total_segments);

    std::deque<u64> window;
    auto process_prime = [&](u64 p) {
        window.push_back(p);
        if (window.size() > 4ULL) {
            window.pop_front();
        }
        if (window.size() != 4ULL) {
            return;
        }
        if ((window[1] - window[0]) != 6ULL ||
            (window[2] - window[1]) != 6ULL ||
            (window[3] - window[2]) != 6ULL) {
            return;
        }

        const u64 n = window[0] + 9ULL;
        if (n > limit) {
            return;
        }
        if (!is_engineers_paradise_candidate(n, factor_primes)) {
            return;
        }

        result.paradises.push_back(n);
        result.sum += n;
        std::cout << "Found paradise #" << result.paradises.size() << ": " << n << '\n';
    };

    if (prime_limit >= 2ULL) {
        process_prime(2ULL);
    }

    while (segment_low <= prime_limit && result.paradises.size() < target_count) {
        std::vector<std::pair<u64, u64>> batch_ranges;
        batch_ranges.reserve(threads);
        for (unsigned t = 0; t < threads && segment_low <= prime_limit; ++t) {
            const u64 low = segment_low;
            const u64 high = std::min<u64>(prime_limit, low + segment_span - 1ULL);
            batch_ranges.emplace_back(low, high);
            segment_low = high + 1ULL;
        }

        std::vector<std::vector<u64>> batch_primes(batch_ranges.size());
        if (batch_ranges.size() == 1ULL) {
            batch_primes[0] = sieve_segment(batch_ranges[0].first, batch_ranges[0].second, base_primes);
        } else {
            std::vector<std::thread> workers;
            workers.reserve(batch_ranges.size());
            for (std::size_t i = 0; i < batch_ranges.size(); ++i) {
                workers.emplace_back([&, i]() {
                    batch_primes[i] =
                        sieve_segment(batch_ranges[i].first, batch_ranges[i].second, base_primes);
                });
            }
            for (std::thread& worker : workers) {
                worker.join();
            }
        }

        for (const std::vector<u64>& primes : batch_primes) {
            for (const u64 p : primes) {
                process_prime(p);
                if (result.paradises.size() >= target_count) {
                    break;
                }
            }
            if (result.paradises.size() >= target_count) {
                break;
            }
        }
    }

    return result;
}

bool check_known_first_two(const SearchResult& result, const Options& options) {
    if (options.limit < kKnownParadise2 || options.target_count < 2U) {
        return true;
    }
    if (result.paradises.size() < 2ULL) {
        std::cerr << "Known-value checkpoint failed: fewer than two paradises found.\n";
        return false;
    }
    if (result.paradises[0] != kKnownParadise1 || result.paradises[1] != kKnownParadise2) {
        std::cerr << "Known-value checkpoint failed.\n";
        std::cerr << "Expected first two paradises: "
                  << kKnownParadise1 << ", " << kKnownParadise2 << '\n';
        std::cerr << "Observed: "
                  << result.paradises[0] << ", " << result.paradises[1] << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const auto start_time = std::chrono::steady_clock::now();

    const u64 sieve_max = options.limit + 9ULL;
    const int root = static_cast<int>(isqrt_u64(sieve_max)) + 1;
    const std::vector<int> base_primes = sieve_primes(root);

    if (options.run_checkpoints) {
        if (!check_practical_criterion(base_primes)) {
            return 1;
        }
        std::cout << "Checkpoint passed: practical criterion matches brute force up to "
                  << kPracticalCrossCheckMax << ".\n";
    }

    SearchResult result = find_paradises(options.limit,
                                         options.target_count,
                                         options.segment_span,
                                         options.allow_multithreading,
                                         options.requested_threads,
                                         base_primes,
                                         base_primes);

    if (options.run_checkpoints) {
        if (!check_known_first_two(result, options)) {
            return 1;
        }
        if (options.limit >= kKnownParadise2 && options.target_count >= 2U) {
            std::cout << "Checkpoint passed: first two paradises match known values.\n";
        }
    }

    if (result.paradises.size() < options.target_count) {
        std::cerr << "Only found " << result.paradises.size()
                  << " paradises up to limit " << options.limit
                  << ". Increase --limit.\n";
        return 1;
    }

    const auto end_time = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "Sum of first " << options.target_count
              << " engineers' paradises: " << result.sum << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " seconds\n";

    return 0;
}
