#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr u64 kDefaultM = 64ULL;
constexpr u64 kDefaultN = 10'000'000'000'000'000ULL;

struct Checkpoint {
    u64 m = 0ULL;
    u64 n = 0ULL;
    u64 expected = 0ULL;
};

constexpr Checkpoint kPublishedCheckpoints[] = {
    {64ULL, 64ULL, 1'263ULL},
    {12ULL, 345ULL, 1'998ULL},
    {32ULL, 1'000'000'000'000'000ULL, 13'826'382'602'124'302ULL},
};

struct Options {
    u64 m = kDefaultM;
    u64 n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
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
        if (parse_u64_after_prefix(arg, "--m=", parsed_u64)) {
            options.m = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            options.n = parsed_u64;
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

    if (options.m == 0ULL) {
        std::cerr << "--m must be positive.\n";
        return false;
    }
    if (options.n == 0ULL) {
        std::cerr << "--n must be positive.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    constexpr std::size_t kMinUnitsForParallel = 8ULL;
    if (!allow_multithreading || workload_units < kMinUnitsForParallel) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

u64 bounded_lcm_or_zero(const u64 a, const u64 b, const u64 limit) {
    const u64 g = std::gcd(a, b);
    const u64 reduced = a / g;
    if (reduced > limit / b) {
        return 0ULL;
    }
    return reduced * b;
}

std::vector<u64> build_minimal_forbidden_divisors(const u64 m, const u64 d, const u64 n) {
    std::vector<u64> all_q;
    all_q.reserve(static_cast<std::size_t>(m - d));
    for (u64 e = d + 1ULL; e <= m; ++e) {
        const u64 q = e / std::gcd(e, d);
        if (q > 1ULL && q <= n) {
            all_q.push_back(q);
        }
    }

    std::sort(all_q.begin(), all_q.end());
    all_q.erase(std::unique(all_q.begin(), all_q.end()), all_q.end());

    std::vector<u64> minimal;
    minimal.reserve(all_q.size());
    for (const u64 q : all_q) {
        bool redundant = false;
        for (const u64 kept : minimal) {
            if (q % kept == 0ULL) {
                redundant = true;
                break;
            }
        }
        if (!redundant) {
            minimal.push_back(q);
        }
    }

    std::sort(minimal.begin(), minimal.end(), std::greater<u64>());
    return minimal;
}

void inclusion_exclusion_dfs(const std::vector<u64>& divisors,
                             const u64 n,
                             const std::size_t start_idx,
                             const u64 current_lcm,
                             const i128 sign,
                             i128& accumulator) {
    for (std::size_t i = start_idx; i < divisors.size(); ++i) {
        const u64 next_lcm = bounded_lcm_or_zero(current_lcm, divisors[i], n);
        if (next_lcm == 0ULL || next_lcm > n) {
            continue;
        }

        accumulator += sign * static_cast<i128>(n / next_lcm);
        inclusion_exclusion_dfs(divisors, n, i + 1ULL, next_lcm, -sign, accumulator);
    }
}

u64 count_not_divisible_by_any(const u64 n, const std::vector<u64>& divisors) {
    // Inclusion-exclusion over forbidden divisibility constraints.
    i128 acc = static_cast<i128>(n);
    inclusion_exclusion_dfs(divisors, n, 0ULL, 1ULL, -1, acc);
    return static_cast<u64>(acc);
}

u64 count_row_contribution(const u64 m, const u64 n, const u64 d) {
    // x = d*b is counted in row d iff no e>d divides x; this turns into q \nmid b tests.
    const std::vector<u64> forbidden = build_minimal_forbidden_divisors(m, d, n);
    return count_not_divisible_by_any(n, forbidden);
}

u64 solve(const u64 m, const u64 n, const bool allow_multithreading, const unsigned requested_threads) {
    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, static_cast<std::size_t>(m));
    if (threads <= 1U) {
        u64 total = 0ULL;
        for (u64 d = 1ULL; d <= m; ++d) {
            total += count_row_contribution(m, n, d);
        }
        return total;
    }

    std::vector<u64> row_counts(static_cast<std::size_t>(m) + 1ULL, 0ULL);
    std::atomic<u64> next_d(1ULL);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned t = 0U; t < threads; ++t) {
        pool.emplace_back([&]() {
            while (true) {
                const u64 d = next_d.fetch_add(1ULL, std::memory_order_relaxed);
                if (d > m) {
                    break;
                }
                row_counts[static_cast<std::size_t>(d)] = count_row_contribution(m, n, d);
            }
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    u64 total = 0ULL;
    for (u64 d = 1ULL; d <= m; ++d) {
        total += row_counts[static_cast<std::size_t>(d)];
    }
    return total;
}

u64 brute_force_distinct_products(const u64 m, const u64 n) {
    const u64 max_product = m * n;
    std::vector<unsigned char> seen(static_cast<std::size_t>(max_product) + 1ULL, 0U);
    for (u64 i = 1ULL; i <= m; ++i) {
        for (u64 j = 1ULL; j <= n; ++j) {
            seen[static_cast<std::size_t>(i * j)] = 1U;
        }
    }

    u64 count = 0ULL;
    for (u64 x = 1ULL; x <= max_product; ++x) {
        count += static_cast<u64>(seen[static_cast<std::size_t>(x)] != 0U);
    }
    return count;
}

bool run_checkpoints() {
    for (const Checkpoint& cp : kPublishedCheckpoints) {
        const u64 got = solve(cp.m, cp.n, false, 1U);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed: P(" << cp.m << ", " << cp.n << ") expected "
                      << cp.expected << " but got " << got << ".\n";
            return false;
        }
    }

    const Checkpoint brute_checks[] = {
        {3ULL, 4ULL, 8ULL},
        {8ULL, 50ULL, 0ULL},
        {16ULL, 80ULL, 0ULL},
    };

    for (const Checkpoint& cp : brute_checks) {
        const u64 got = solve(cp.m, cp.n, false, 1U);
        const u64 brute = brute_force_distinct_products(cp.m, cp.n);
        if (got != brute) {
            std::cerr << "Brute-force cross-check failed: P(" << cp.m << ", " << cp.n
                      << ") optimized=" << got << " brute=" << brute << ".\n";
            return false;
        }
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
        return 1;
    }

    const u64 answer =
        solve(options.m, options.n, options.allow_multithreading, options.requested_threads);
    std::cout << answer << '\n';

    return 0;
}
