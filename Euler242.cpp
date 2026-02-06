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
using u128 = __uint128_t;

constexpr u64 kDefaultLimit = 1'000'000'000'000ULL;

struct Options {
    u64 limit = kDefaultLimit;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string result;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        result.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }

    std::reverse(result.begin(), result.end());
    return result;
}

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
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
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

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload) {
    if (!allow_multithreading || workload < 8ULL) {
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

bool binom_is_odd(const u64 n, const u64 k) {
    return k <= n && ((k & ~n) == 0ULL);
}

u128 binom_exact(const u64 n, u64 k) {
    if (k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }

    u128 value = 1;
    for (u64 i = 1; i <= k; ++i) {
        value = (value * static_cast<u128>(n - k + i)) / static_cast<u128>(i);
    }
    return value;
}

u128 f_exact(const u64 n, const u64 k) {
    const u64 odd_count = (n + 1ULL) / 2ULL;
    const u64 even_count = n / 2ULL;

    u128 total = 0;
    for (u64 odd_taken = 1ULL; odd_taken <= k && odd_taken <= odd_count; odd_taken += 2ULL) {
        const u64 even_taken = k - odd_taken;
        if (even_taken > even_count) {
            continue;
        }
        total += binom_exact(odd_count, odd_taken) * binom_exact(even_count, even_taken);
    }
    return total;
}

bool f_is_odd_direct(const u64 n, const u64 k) {
    if (k > n) {
        return false;
    }

    const u64 odd_count = (n + 1ULL) / 2ULL;
    const u64 even_count = n / 2ULL;

    bool parity = false;
    for (u64 odd_taken = 1ULL; odd_taken <= k && odd_taken <= odd_count; odd_taken += 2ULL) {
        const u64 even_taken = k - odd_taken;
        if (even_taken > even_count) {
            continue;
        }

        if (binom_is_odd(odd_count, odd_taken) && binom_is_odd(even_count, even_taken)) {
            parity = !parity;
        }
    }
    return parity;
}

bool f_is_odd_formula(const u64 n, const u64 k) {
    if ((n & 1ULL) == 0ULL || (k & 1ULL) == 0ULL || k > n) {
        return false;
    }

    const u64 m = (n - 1ULL) / 2ULL;
    if ((m & 1ULL) != 0ULL) {
        return false;
    }

    const u64 r = (k - 1ULL) / 2ULL;
    return binom_is_odd(m, r);
}

u64 count_odd_triplets_bruteforce(const u64 limit) {
    u64 total = 0;
    for (u64 n = 1ULL; n <= limit; n += 2ULL) {
        for (u64 k = 1ULL; k <= n; k += 2ULL) {
            if (f_is_odd_direct(n, k)) {
                ++total;
            }
        }
    }
    return total;
}

const std::vector<u128>& powers_of_three() {
    static const std::vector<u128> table = [] {
        std::vector<u128> values(65, 0);
        values[0] = 1;
        for (std::size_t i = 1; i < values.size(); ++i) {
            values[i] = values[i - 1] * static_cast<u128>(3);
        }
        return values;
    }();
    return table;
}

u128 weighted_popcount_prefix_sum(const u64 n) {
    const auto& pow3 = powers_of_three();

    u128 result = 0;
    unsigned ones_so_far = 0U;

    for (int bit = 63; bit >= 0; --bit) {
        if (((n >> bit) & 1ULL) == 0ULL) {
            continue;
        }

        result += (static_cast<u128>(1) << ones_so_far) * pow3[static_cast<std::size_t>(bit)];
        ++ones_so_far;
    }

    result += (static_cast<u128>(1) << ones_so_far);
    return result;
}

u128 solve(const u64 limit) {
    if (limit == 0ULL) {
        return 0;
    }

    // n = 2m+1 and only even m contribute, so m = 2t.
    const u64 m_max = (limit - 1ULL) / 2ULL;
    const u64 t_max = m_max / 2ULL;
    return weighted_popcount_prefix_sum(t_max);
}

bool run_checkpoints(const Options& options) {
    bool ok = true;

    if (f_exact(5, 3) != static_cast<u128>(4)) {
        std::cerr << "Checkpoint failed: f(5,3) should be 4.\n";
        ok = false;
    }

    const u64 sample_count = count_odd_triplets_bruteforce(10ULL);
    if (sample_count != 5ULL) {
        std::cerr << "Checkpoint failed: odd-triplets up to n=10 should be 5.\n";
        ok = false;
    }

    const u64 parity_check_max_n = 1001ULL;
    const std::size_t workload = static_cast<std::size_t>((parity_check_max_n + 1ULL) / 2ULL);
    const unsigned thread_count =
        choose_thread_count(options.allow_multithreading, options.requested_threads, workload);

    std::atomic<bool> parity_ok(true);
    std::atomic<u64> next_n(1ULL);

    auto worker = [&]() {
        while (parity_ok.load(std::memory_order_relaxed)) {
            const u64 n = next_n.fetch_add(2ULL, std::memory_order_relaxed);
            if (n > parity_check_max_n) {
                break;
            }

            for (u64 k = 1ULL; k <= n; k += 2ULL) {
                if (f_is_odd_direct(n, k) != f_is_odd_formula(n, k)) {
                    parity_ok.store(false, std::memory_order_relaxed);
                    return;
                }
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(thread_count);
    for (unsigned i = 0U; i < thread_count; ++i) {
        pool.emplace_back(worker);
    }
    for (auto& th : pool) {
        th.join();
    }

    if (!parity_ok.load()) {
        std::cerr << "Checkpoint failed: parity identity mismatch.\n";
        ok = false;
    }

    u128 running = 0;
    constexpr u64 kPrefixCheckMax = 200000ULL;
    for (u64 x = 0ULL; x <= kPrefixCheckMax; ++x) {
        running += (static_cast<u128>(1) << __builtin_popcountll(x));
        if (weighted_popcount_prefix_sum(x) != running) {
            std::cerr << "Checkpoint failed: weighted prefix mismatch at x=" << x << ".\n";
            ok = false;
            break;
        }
    }

    if (ok) {
        std::cout << "Checkpoints passed (threads=" << thread_count << ").\n";
    }

    return ok;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const u128 answer = solve(options.limit);
    std::cout << "Answer: " << to_string_u128(answer) << '\n';
    return 0;
}
