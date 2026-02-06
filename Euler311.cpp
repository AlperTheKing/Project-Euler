#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u16 = std::uint16_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultLimit = 10'000'000'000ULL;
constexpr u64 kDefaultBlockSpan = 4'000'000ULL;
constexpr u64 kCheckpointLimit1 = 10'000ULL;
constexpr u64 kCheckpointExpected1 = 49ULL;
constexpr u64 kCheckpointLimit2 = 1'000'000ULL;
constexpr u64 kCheckpointExpected2 = 38'239ULL;
constexpr u64 kThreadConsistencyLimit = 4'000'000ULL;

struct Options {
    u64 limit = kDefaultLimit;
    u64 block_span = kDefaultBlockSpan;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
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
        if (parse_u64_after_prefix(arg, "--block=", parsed_u64)) {
            options.block_span = parsed_u64;
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

    if (options.block_span == 0ULL) {
        std::cerr << "--block must be at least 1.\n";
        return false;
    }

    return true;
}

u64 isqrt_u64(u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

u64 isqrt_ceil_u64(u64 n) {
    const u64 r = isqrt_u64(n);
    return (r * r == n) ? r : (r + 1ULL);
}

bool is_twice_square(u64 n) {
    if ((n & 1ULL) != 0ULL) {
        return false;
    }
    const u64 m = n >> 1ULL;
    const u64 r = isqrt_u64(m);
    return r * r == m;
}

u64 n_choose_2(u64 n) {
    return (n < 2ULL) ? 0ULL : (n * (n - 1ULL) / 2ULL);
}

u64 n_choose_3(u64 n) {
    if (n < 3ULL) {
        return 0ULL;
    }
    const u128 numerator = static_cast<u128>(n) * (n - 1ULL) * (n - 2ULL);
    return static_cast<u64>(numerator / 6ULL);
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

u64 process_block(u64 low,
                  u64 high,
                  std::vector<u16>& counts,
                  std::vector<u32>& touched) {
    touched.clear();

    for (u64 y = 1ULL;; ++y) {
        const u64 y2 = y * y;
        const u64 min_n = y2 + (y + 1ULL) * (y + 1ULL);
        if (min_n > high) {
            break;
        }

        u64 x_min = y + 1ULL;
        const u64 x_min_n = y2 + x_min * x_min;
        if (x_min_n < low) {
            x_min = isqrt_ceil_u64(low - y2);
            if (x_min <= y) {
                x_min = y + 1ULL;
            }
        }

        const u64 x_max = isqrt_u64(high - y2);
        if (x_max < x_min) {
            continue;
        }

        u64 x2 = x_min * x_min;
        for (u64 x = x_min; x <= x_max; ++x) {
            const u64 n = y2 + x2;
            const u32 idx = static_cast<u32>(n - low);
            u16& count = counts[static_cast<std::size_t>(idx)];
            if (count == 0U) {
                touched.push_back(idx);
            }
            if (count == std::numeric_limits<u16>::max()) {
                throw std::runtime_error("Representation count overflowed uint16_t.");
            }
            ++count;
            x2 += (2ULL * x + 1ULL);
        }
    }

    u64 block_sum = 0ULL;
    for (const u32 idx : touched) {
        u16& count = counts[static_cast<std::size_t>(idx)];
        const u64 m = static_cast<u64>(count);
        if (m >= 2ULL) {
            const u64 n = low + static_cast<u64>(idx);
            block_sum += n_choose_3(m);
            if (is_twice_square(n)) {
                block_sum += n_choose_2(m);
            }
        }
        count = 0U;
    }

    return block_sum;
}

u64 solve_biclinic(u64 limit,
                   bool allow_multithreading,
                   unsigned requested_threads,
                   u64 block_span) {
    if (limit < 4ULL) {
        return 0ULL;
    }

    const u64 max_t = limit / 4ULL;
    const u64 block_count_u64 = (max_t + block_span - 1ULL) / block_span;
    if (block_count_u64 > static_cast<u64>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("Too many blocks for this platform.");
    }
    const std::size_t block_count = static_cast<std::size_t>(block_count_u64);
    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, block_count);

    std::vector<u64> partial(threads, 0ULL);
    auto worker = [&](unsigned tid) {
        std::vector<u16> counts(static_cast<std::size_t>(block_span), 0U);
        std::vector<u32> touched;
        touched.reserve(static_cast<std::size_t>(block_span / 4ULL + 1024ULL));

        u64 local_sum = 0ULL;
        for (std::size_t bi = tid; bi < block_count; bi += threads) {
            const u64 low = 1ULL + static_cast<u64>(bi) * block_span;
            const u64 high = std::min<u64>(max_t, low + block_span - 1ULL);
            local_sum += process_block(low, high, counts, touched);
        }
        partial[tid] = local_sum;
    };

    std::vector<std::thread> pool;
    pool.reserve(threads > 0U ? threads - 1U : 0U);
    for (unsigned t = 1U; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    worker(0U);
    for (std::thread& thread : pool) {
        thread.join();
    }

    u64 total = 0ULL;
    for (const u64 x : partial) {
        total += x;
    }
    return total;
}

bool run_checkpoints(const Options& options) {
    const u64 sample1 = solve_biclinic(
        kCheckpointLimit1, options.allow_multithreading, options.requested_threads, options.block_span);
    if (sample1 != kCheckpointExpected1) {
        std::cerr << "Checkpoint failed for N=" << kCheckpointLimit1
                  << ": expected " << kCheckpointExpected1
                  << ", got " << sample1 << '\n';
        return false;
    }
    std::cout << "Checkpoint passed: B(" << kCheckpointLimit1 << ") = " << sample1 << '\n';

    const u64 sample2 = solve_biclinic(
        kCheckpointLimit2, options.allow_multithreading, options.requested_threads, options.block_span);
    if (sample2 != kCheckpointExpected2) {
        std::cerr << "Checkpoint failed for N=" << kCheckpointLimit2
                  << ": expected " << kCheckpointExpected2
                  << ", got " << sample2 << '\n';
        return false;
    }
    std::cout << "Checkpoint passed: B(" << kCheckpointLimit2 << ") = " << sample2 << '\n';

    if (options.allow_multithreading) {
        const u64 multi = solve_biclinic(
            kThreadConsistencyLimit, true, options.requested_threads, options.block_span);
        const u64 single = solve_biclinic(
            kThreadConsistencyLimit, false, options.requested_threads, options.block_span);
        if (multi != single) {
            std::cerr << "Thread consistency failed at N=" << kThreadConsistencyLimit
                      << ": multi-thread=" << multi
                      << ", single-thread=" << single << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: thread consistency at N="
                  << kThreadConsistencyLimit << '\n';
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

    try {
        const auto start_time = std::chrono::steady_clock::now();

        if (options.run_checkpoints && !run_checkpoints(options)) {
            return 1;
        }

        const u64 result = solve_biclinic(options.limit,
                                          options.allow_multithreading,
                                          options.requested_threads,
                                          options.block_span);

        const auto end_time = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end_time - start_time;

        std::cout << "B(" << options.limit << ") = " << result << '\n';
        std::cout << "Elapsed: " << elapsed.count() << " seconds\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}

