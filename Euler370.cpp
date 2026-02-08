#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

constexpr u64 kDefaultLimit = 25'000'000'000'000ULL;
constexpr u64 kThreadConsistencyLimit = 100'000'000'000ULL;

struct Checkpoint {
    u64 limit;
    u64 expected;
    const char* label;
};

constexpr std::array<Checkpoint, 5> kStatementCheckpoints = {{
    {100ULL, 42ULL, "G(100)"},
    {1'000ULL, 532ULL, "G(1000)"},
    {10'000ULL, 6'427ULL, "G(10^4)"},
    {100'000ULL, 75'243ULL, "G(10^5)"},
    {1'000'000ULL, 861'805ULL, "G(10^6)"},
}};

struct Options {
    u64 limit = kDefaultLimit;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

u64 isqrt_u128(const u128 x) {
    long double root = std::sqrt(static_cast<long double>(x));
    u64 value = static_cast<u64>(root);

    while (static_cast<u128>(value + 1ULL) * static_cast<u128>(value + 1ULL) <= x) {
        ++value;
    }
    while (static_cast<u128>(value) * static_cast<u128>(value) > x) {
        --value;
    }
    return value;
}

u64 isqrt_u64(const u64 x) {
    return isqrt_u128(static_cast<u128>(x));
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }

    std::string digits;
    while (value > 0U) {
        const int digit = static_cast<int>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
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

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--limit=", parsed_u64)) {
            options.limit = parsed_u64;
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.limit < 3ULL) {
        std::cerr << "--limit must be at least 3.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const u64 work_items) {
    if (!allow_multithreading || work_items <= 1ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    threads = std::min<unsigned>(threads, static_cast<unsigned>(work_items));
    return std::max(1U, threads);
}

class GeometricTriangleCounter {
public:
    explicit GeometricTriangleCounter(const u64 max_limit)
        : p_max_capacity_(isqrt_u64(max_limit / 3ULL)) {
        build_smallest_prime_factors(static_cast<int>(p_max_capacity_));
    }

    u128 count(const u64 limit, const unsigned requested_threads) const {
        const u64 p_max = isqrt_u64(limit / 3ULL);
        if (p_max == 0ULL) {
            return 0U;
        }
        if (p_max > p_max_capacity_) {
            return 0U;
        }

        const unsigned threads =
            choose_thread_count(true, requested_threads, p_max);

        if (threads == 1U) {
            return count_range(limit, 1ULL, p_max);
        }

        std::atomic<u64> next_p{1ULL};
        std::vector<u128> partial(threads, 0U);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0U; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                u128 local = 0U;
                while (true) {
                    const u64 p = next_p.fetch_add(1ULL, std::memory_order_relaxed);
                    if (p > p_max) {
                        break;
                    }
                    local += count_for_single_p(limit, p);
                }
                partial[tid] = local;
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        u128 total = 0U;
        for (const u128 part : partial) {
            total += part;
        }
        return total;
    }

    static u64 max_q_triangle(const u64 p) {
        return q_max_triangle(p);
    }

    static u64 max_q_for_norm(const u64 p, const u64 limit) {
        return q_max_for_norm(p, limit);
    }

private:
    std::vector<int> smallest_prime_factor_;
    u64 p_max_capacity_ = 0ULL;

    static bool triangle_ok(const u64 p, const u64 q) {
        return static_cast<u128>(q) * static_cast<u128>(q) <
               static_cast<u128>(p) * static_cast<u128>(p) +
                   static_cast<u128>(p) * static_cast<u128>(q);
    }

    static bool norm_ok(const u64 p, const u64 q, const u64 limit) {
        return static_cast<u128>(q) * static_cast<u128>(q) +
                   static_cast<u128>(p) * static_cast<u128>(q) +
                   static_cast<u128>(p) * static_cast<u128>(p) <=
               static_cast<u128>(limit);
    }

    static long double phi_value() {
        static const long double kPhi =
            (1.0L + std::sqrt(5.0L)) / 2.0L;
        return kPhi;
    }

    static u64 q_max_triangle(const u64 p) {
        u64 q = static_cast<u64>(phi_value() * static_cast<long double>(p));
        while (q > 0ULL && !triangle_ok(p, q)) {
            --q;
        }
        while (triangle_ok(p, q + 1ULL)) {
            ++q;
        }
        return q;
    }

    static u64 q_max_for_norm(const u64 p, const u64 limit) {
        if (static_cast<u128>(3U) * static_cast<u128>(p) * static_cast<u128>(p) >
            static_cast<u128>(limit)) {
            return 0ULL;
        }

        const u128 discriminant = static_cast<u128>(4U) * static_cast<u128>(limit) -
                                  static_cast<u128>(3U) * static_cast<u128>(p) *
                                      static_cast<u128>(p);
        const u64 root = isqrt_u128(discriminant);

        i64 q = (static_cast<i64>(root) - static_cast<i64>(p)) / 2LL;
        if (q < 0LL) {
            q = 0LL;
        }

        while (q > 0LL && !norm_ok(p, static_cast<u64>(q), limit)) {
            --q;
        }
        while (norm_ok(p, static_cast<u64>(q) + 1ULL, limit)) {
            ++q;
        }
        return static_cast<u64>(q);
    }

    void build_smallest_prime_factors(const int limit) {
        smallest_prime_factor_.assign(static_cast<std::size_t>(limit) + 1ULL, 0);
        if (limit >= 1) {
            smallest_prime_factor_[1] = 1;
        }

        for (int i = 2; i <= limit; ++i) {
            if (smallest_prime_factor_[static_cast<std::size_t>(i)] == 0) {
                smallest_prime_factor_[static_cast<std::size_t>(i)] = i;
                if (static_cast<i64>(i) * static_cast<i64>(i) <= limit) {
                    for (int j = i * i; j <= limit; j += i) {
                        if (smallest_prime_factor_[static_cast<std::size_t>(j)] == 0) {
                            smallest_prime_factor_[static_cast<std::size_t>(j)] = i;
                        }
                    }
                }
            }
        }
    }

    int extract_distinct_primes(u64 value, std::array<u32, 8>& primes) const {
        int count = 0;
        while (value > 1ULL) {
            const u32 prime = static_cast<u32>(
                smallest_prime_factor_[static_cast<std::size_t>(value)]);
            primes[static_cast<std::size_t>(count)] = prime;
            ++count;
            while (value % static_cast<u64>(prime) == 0ULL) {
                value /= static_cast<u64>(prime);
            }
        }
        return count;
    }

    static int build_squarefree_divisors_mu(const std::array<u32, 8>& primes,
                                            const int prime_count,
                                            std::array<u32, 128>& divisors,
                                            std::array<int8_t, 128>& mus) {
        divisors[0] = 1U;
        mus[0] = 1;
        int count = 1;

        for (int i = 0; i < prime_count; ++i) {
            const int current = count;
            const u32 p = primes[static_cast<std::size_t>(i)];
            for (int j = 0; j < current; ++j) {
                divisors[static_cast<std::size_t>(count)] =
                    divisors[static_cast<std::size_t>(j)] * p;
                mus[static_cast<std::size_t>(count)] =
                    static_cast<int8_t>(-mus[static_cast<std::size_t>(j)]);
                ++count;
            }
        }

        return count;
    }

    static u64 count_coprime_in_interval(
        const u64 p,
        const u64 left,
        const u64 right,
        const std::array<u32, 128>& divisors,
        const std::array<int8_t, 128>& mus,
        const int divisor_count) {
        if (left > right) {
            return 0ULL;
        }

        const u64 span = right - left + 1ULL;
        if (span <= 8ULL) {
            u64 count = 0ULL;
            for (u64 q = left; q <= right; ++q) {
                if (std::gcd(p, q) == 1ULL) {
                    ++count;
                }
            }
            return count;
        }

        i64 total = 0LL;
        for (int i = 0; i < divisor_count; ++i) {
            const u64 d = static_cast<u64>(divisors[static_cast<std::size_t>(i)]);
            total += static_cast<i64>(mus[static_cast<std::size_t>(i)]) *
                     static_cast<i64>(right / d - (left - 1ULL) / d);
        }
        return static_cast<u64>(total);
    }

    u128 count_for_single_p(const u64 limit, const u64 p) const {
        const u64 q_start = p;
        const u64 q_end =
            std::min(q_max_triangle(p), q_max_for_norm(p, limit));
        if (q_end < q_start) {
            return 0U;
        }

        std::array<u32, 8> primes{};
        std::array<u32, 128> divisors{};
        std::array<int8_t, 128> mus{};

        const int prime_count = extract_distinct_primes(p, primes);
        const int divisor_count =
            build_squarefree_divisors_mu(primes, prime_count, divisors, mus);

        u128 total = 0U;
        u64 left = q_start;

        while (left <= q_end) {
            const u64 norm_left = static_cast<u64>(
                static_cast<u128>(left) * static_cast<u128>(left) +
                static_cast<u128>(p) * static_cast<u128>(left) +
                static_cast<u128>(p) * static_cast<u128>(p));
            const u64 quotient = limit / norm_left;
            if (quotient == 0ULL) {
                break;
            }

            const u64 norm_limit = limit / quotient;
            const u64 right =
                std::min(q_end, q_max_for_norm(p, norm_limit));

            const u64 coprime_count = count_coprime_in_interval(
                p, left, right, divisors, mus, divisor_count);
            total += static_cast<u128>(coprime_count) *
                     static_cast<u128>(quotient);

            left = right + 1ULL;
        }

        return total;
    }

    u128 count_range(const u64 limit, const u64 start_p, const u64 end_p) const {
        u128 total = 0U;
        for (u64 p = start_p; p <= end_p; ++p) {
            total += count_for_single_p(limit, p);
        }
        return total;
    }
};

u64 brute_force_count(const u64 limit) {
    const u64 p_max = isqrt_u64(limit / 3ULL);
    u64 total = 0ULL;

    for (u64 p = 1ULL; p <= p_max; ++p) {
        const u64 q_end =
            std::min(GeometricTriangleCounter::max_q_triangle(p),
                     GeometricTriangleCounter::max_q_for_norm(p, limit));
        if (q_end < p) {
            continue;
        }

        for (u64 q = p; q <= q_end; ++q) {
            if (std::gcd(p, q) != 1ULL) {
                continue;
            }

            const u64 norm = static_cast<u64>(
                static_cast<u128>(p) * static_cast<u128>(p) +
                static_cast<u128>(p) * static_cast<u128>(q) +
                static_cast<u128>(q) * static_cast<u128>(q));
            total += limit / norm;
        }
    }

    return total;
}

bool check_equal_u128(const char* label, const u128 got, const u128 expected) {
    if (got != expected) {
        std::cerr << "Checkpoint failed for " << label << ": expected "
                  << to_string_u128(expected) << ", got " << to_string_u128(got)
                  << '\n';
        return false;
    }
    std::cout << "Checkpoint OK: " << label << " = " << to_string_u128(got)
              << '\n';
    return true;
}

bool run_checkpoints(const GeometricTriangleCounter& counter,
                     const unsigned threads) {
    for (const Checkpoint checkpoint : kStatementCheckpoints) {
        const u128 got = counter.count(checkpoint.limit, 1U);
        if (!check_equal_u128(checkpoint.label, got,
                              static_cast<u128>(checkpoint.expected))) {
            return false;
        }
    }

    constexpr std::array<u64, 6> kBruteLimits = {
        3ULL, 17ULL, 111ULL, 999ULL, 12'345ULL, 54'321ULL};
    for (const u64 limit : kBruteLimits) {
        const u128 fast = counter.count(limit, 1U);
        const u128 brute = static_cast<u128>(brute_force_count(limit));
        const std::string label = "fast-vs-brute G(" + std::to_string(limit) + ")";
        if (!check_equal_u128(label.c_str(), fast, brute)) {
            return false;
        }
    }

    if (threads > 1U) {
        const u128 single = counter.count(kThreadConsistencyLimit, 1U);
        const u128 multi = counter.count(kThreadConsistencyLimit, threads);
        if (!check_equal_u128("thread consistency",
                              multi,
                              single)) {
            return false;
        }
    }

    std::cout << "Validation checkpoints passed.\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    u64 max_limit_for_solver = options.limit;
    if (options.run_checkpoints) {
        for (const Checkpoint checkpoint : kStatementCheckpoints) {
            max_limit_for_solver = std::max(max_limit_for_solver, checkpoint.limit);
        }
        max_limit_for_solver =
            std::max(max_limit_for_solver, kThreadConsistencyLimit);
    }

    const GeometricTriangleCounter counter(max_limit_for_solver);
    const u64 p_max = isqrt_u64(options.limit / 3ULL);
    const unsigned threads =
        choose_thread_count(options.allow_multithreading,
                            options.requested_threads,
                            p_max);

    if (options.run_checkpoints) {
        if (!run_checkpoints(counter, threads)) {
            return 1;
        }
    }

    const u128 answer = counter.count(options.limit, threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
