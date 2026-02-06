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
using i128 = __int128_t;

constexpr u64 kDefaultLimitExclusive = 10'000'000'000'000'000ULL;
constexpr u64 kSampleLimitExclusive = 10'000'000ULL;
constexpr u64 kSampleExpected = 2'327'192ULL;
constexpr u64 kBruteCheckpointLimitExclusive = 200'003ULL;

struct Options {
    u64 limit_exclusive = kDefaultLimitExclusive;
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
            options.limit_exclusive = limit;
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

    return true;
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

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2'000'000ULL) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    return std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

struct SieveData {
    std::vector<std::int8_t> mobius_odd;
    std::vector<int> primes_3mod4;
};

SieveData build_sieve_data(u64 max_u) {
    const std::size_t odd_count = static_cast<std::size_t>(max_u / 2ULL + 1ULL);
    std::vector<std::uint8_t> composite(odd_count, 0U);
    std::vector<int> odd_primes;
    if (max_u >= 3ULL) {
        odd_primes.reserve(static_cast<std::size_t>(max_u / std::max<u64>(
                                                                1ULL,
                                                                static_cast<u64>(
                                                                    std::log(static_cast<long double>(max_u))))));
    }

    for (std::size_t idx = 1; idx < odd_count; ++idx) {
        if (composite[idx] != 0U) {
            continue;
        }

        const u64 p = 2ULL * static_cast<u64>(idx) + 1ULL;
        odd_primes.push_back(static_cast<int>(p));

        if (p > max_u / p) {
            continue;
        }

        std::size_t mark = static_cast<std::size_t>((p * p - 1ULL) / 2ULL);
        while (mark < odd_count) {
            composite[mark] = 1U;
            mark += static_cast<std::size_t>(p);
        }
    }

    SieveData out;
    out.mobius_odd.assign(odd_count, static_cast<std::int8_t>(1));
    out.primes_3mod4.reserve(odd_primes.size() / 2ULL + 1ULL);

    for (const int p_int : odd_primes) {
        const u64 p = static_cast<u64>(p_int);
        if ((p & 3ULL) == 3ULL) {
            out.primes_3mod4.push_back(p_int);
        }

        std::size_t idx = static_cast<std::size_t>((p - 1ULL) / 2ULL);
        while (idx < odd_count) {
            out.mobius_odd[idx] = static_cast<std::int8_t>(-out.mobius_odd[idx]);
            idx += static_cast<std::size_t>(p);
        }

        if (p > max_u / p) {
            continue;
        }
        const u64 p2 = p * p;
        idx = static_cast<std::size_t>((p2 - 1ULL) / 2ULL);
        while (idx < odd_count) {
            out.mobius_odd[idx] = static_cast<std::int8_t>(0);
            idx += static_cast<std::size_t>(p2);
        }
    }

    return out;
}

std::vector<std::int16_t> build_coefficients(const std::vector<std::int8_t>& mobius_odd,
                                             const std::vector<int>& primes_3mod4) {
    std::vector<std::int16_t> coeff(mobius_odd.size(), 0);
    for (std::size_t i = 0; i < mobius_odd.size(); ++i) {
        coeff[i] = static_cast<std::int16_t>(mobius_odd[i]);
    }

    // c(u) = mu(u) + sum_{q|u, q==3(mod4) prime} mu(u/q), for odd u.
    for (const int q_int : primes_3mod4) {
        const std::size_t q = static_cast<std::size_t>(q_int);
        std::size_t n_idx = (q - 1ULL) / 2ULL;
        std::size_t m_idx = 0;
        while (n_idx < coeff.size()) {
            coeff[n_idx] = static_cast<std::int16_t>(coeff[n_idx] + mobius_odd[m_idx]);
            ++m_idx;
            n_idx += q;
        }
    }

    return coeff;
}

i128 accumulate_range(const std::vector<std::int16_t>& coeff,
                     u64 inclusive_limit,
                     std::size_t begin,
                     std::size_t end) {
    i128 local = 0;
    for (std::size_t idx = begin; idx < end; ++idx) {
        const int c = static_cast<int>(coeff[idx]);
        if (c == 0) {
            continue;
        }
        const u64 u = 2ULL * static_cast<u64>(idx) + 1ULL;
        const u64 u2 = u * u;
        // Number of integers <= inclusive_limit / u^2 that are 1 (mod 4).
        const u64 terms = (inclusive_limit / u2 + 3ULL) / 4ULL;
        local += static_cast<i128>(c) * static_cast<i128>(terms);
    }
    return local;
}

u64 count_squarefree_hilbert_upto(u64 inclusive_limit,
                                  bool allow_multithreading,
                                  unsigned requested_threads) {
    const u64 max_u = isqrt_u64(inclusive_limit);
    const SieveData sieve = build_sieve_data(max_u);
    const std::vector<std::int16_t> coeff =
        build_coefficients(sieve.mobius_odd, sieve.primes_3mod4);

    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, coeff.size());
    if (threads == 1) {
        const i128 total = accumulate_range(coeff, inclusive_limit, 0, coeff.size());
        return static_cast<u64>(total);
    }

    std::vector<std::thread> pool;
    std::vector<i128> partial(threads, 0);
    pool.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        const std::size_t begin = coeff.size() * t / threads;
        const std::size_t end = coeff.size() * (t + 1ULL) / threads;
        pool.emplace_back([&, t, begin, end]() {
            partial[t] = accumulate_range(coeff, inclusive_limit, begin, end);
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    i128 total = 0;
    for (const i128 x : partial) {
        total += x;
    }

    return static_cast<u64>(total);
}

u64 count_squarefree_hilbert_below(u64 exclusive_limit,
                                   bool allow_multithreading,
                                   unsigned requested_threads) {
    if (exclusive_limit == 0ULL) {
        return 0ULL;
    }
    return count_squarefree_hilbert_upto(exclusive_limit - 1ULL,
                                         allow_multithreading,
                                         requested_threads);
}

u64 brute_force_count_below(u64 exclusive_limit) {
    if (exclusive_limit == 0ULL) {
        return 0ULL;
    }

    const u64 inclusive_limit = exclusive_limit - 1ULL;
    const u64 root = isqrt_u64(inclusive_limit);

    std::vector<u64> hilbert_squares;
    for (u64 h = 5ULL; h <= root; h += 4ULL) {
        hilbert_squares.push_back(h * h);
    }

    u64 count = 0;
    for (u64 n = 1ULL; n <= inclusive_limit; n += 4ULL) {
        bool ok = true;
        for (const u64 sq : hilbert_squares) {
            if (sq > n) {
                break;
            }
            if (n % sq == 0ULL) {
                ok = false;
                break;
            }
        }
        if (ok) {
            ++count;
        }
    }

    return count;
}

bool run_checkpoints(const Options& options) {
    struct Point {
        u64 limit_exclusive;
        u64 expected;
        const char* label;
    };

    const std::vector<Point> fixed = {
        {10ULL, 3ULL, "count(<10)"},
        {13ULL, 3ULL, "count(<13)"},
        {100ULL, 23ULL, "count(<100)"},
        {1'000ULL, 232ULL, "count(<1000)"},
        {kSampleLimitExclusive, kSampleExpected, "count(<10^7)"}};

    for (const Point& cp : fixed) {
        const u64 got =
            count_squarefree_hilbert_below(cp.limit_exclusive, false, 1U);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed: " << cp.label << ", expected " << cp.expected
                      << ", got " << got << '\n';
            return false;
        }
        std::cout << "Checkpoint OK: " << cp.label << " = " << got << '\n';
    }

    const u64 brute_expected = brute_force_count_below(kBruteCheckpointLimitExclusive);
    const u64 brute_fast = count_squarefree_hilbert_below(
        kBruteCheckpointLimitExclusive,
        options.allow_multithreading,
        options.requested_threads);
    if (brute_expected != brute_fast) {
        std::cerr << "Brute checkpoint failed at <" << kBruteCheckpointLimitExclusive
                  << ": expected " << brute_expected << ", got " << brute_fast << '\n';
        return false;
    }
    std::cout << "Checkpoint OK: brute cross-check <" << kBruteCheckpointLimitExclusive
              << " = " << brute_fast << '\n';

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const auto t0 = std::chrono::steady_clock::now();

    if (options.run_checkpoints) {
        if (!run_checkpoints(options)) {
            return 1;
        }
    }

    const u64 answer = count_squarefree_hilbert_below(options.limit_exclusive,
                                                      options.allow_multithreading,
                                                      options.requested_threads);

    const auto t1 = std::chrono::steady_clock::now();
    const std::chrono::duration<long double> elapsed = t1 - t0;

    std::cout << "Answer: " << answer << '\n';
    std::cout << "count(<" << options.limit_exclusive << ") = " << answer << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " s\n";

    return 0;
}
