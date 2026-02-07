#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using i32 = std::int32_t;
using i64 = std::int64_t;
using i8 = std::int8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = unsigned __int128;

constexpr u32 kDefaultN = 5'000'000;

struct Options {
    u32 n = kDefaultN;
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
    bool run_validation = true;
    u32 frontier_split_x = 128;
    u32 chunk_divisor = 256;
};

struct FloorPair {
    i128 sum1 = 0;
    i128 sum2 = 0;
};

struct IntervalTask {
    u32 q = 0;
    i32 mu_sum = 0;
};

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 50'000) {
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

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 u32& value_out) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::uint64_t parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value_out = static_cast<u32>(parsed);
    return true;
}

bool parse_options(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        if (arg == "--no-validate") {
            options.run_validation = false;
            continue;
        }

        u32 parsed = 0;
        if (parse_unsigned_after_prefix(arg, "--n=", parsed)) {
            if (parsed == 0) {
                std::cerr << "Invalid --n value: must be >= 1\n";
                return false;
            }
            options.n = parsed;
            continue;
        }

        if (parse_unsigned_after_prefix(arg, "--threads=", parsed)) {
            options.requested_threads = static_cast<unsigned>(parsed);
            continue;
        }

        if (parse_unsigned_after_prefix(arg, "--split=", parsed)) {
            if (parsed == 0) {
                std::cerr << "Invalid --split value: must be >= 1\n";
                return false;
            }
            options.frontier_split_x = parsed;
            continue;
        }

        if (parse_unsigned_after_prefix(arg, "--chunk-div=", parsed)) {
            if (parsed == 0) {
                std::cerr << "Invalid --chunk-div value: must be >= 1\n";
                return false;
            }
            options.chunk_divisor = parsed;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

inline i64 floor_div(i64 a, i64 b) {
    return (a >= 0) ? (a / b) : -(((-a) + b - 1) / b);
}

u32 euclid_steps(u32 x, u32 y) {
    u32 steps = 0;
    while (y != 0) {
        const u32 r = x % y;
        x = y;
        y = r;
        ++steps;
    }
    return steps;
}

u64 brute_force_s(u32 n) {
    u64 total = 0;
    for (u32 x = 1; x <= n; ++x) {
        for (u32 y = 1; y <= n; ++y) {
            total += euclid_steps(x, y);
        }
    }
    return total;
}

FloorPair floor_sum_pair(i64 p, i64 q, i64 m, i64 s1, i64 s2) {
    FloorPair result{};
    i64 sign = 1;

    while (true) {
        i64 t = floor_div(m, q);
        m -= t * q;
        result.sum1 += static_cast<i128>(sign) * t * s1;
        result.sum2 += static_cast<i128>(sign) * t * s2;

        t = floor_div(p, q);
        p -= t * q;
        result.sum1 += static_cast<i128>(sign) * t * s1 * (s1 + 1) / 2;
        result.sum2 += static_cast<i128>(sign) * t * s2 * (s2 + 1) / 2;

        if (p == 0) {
            return result;
        }

        t = (p * s1 + m) / q;
        result.sum1 += static_cast<i128>(sign) * s1 * t;
        s1 = t;

        t = (p * s2 + m) / q;
        result.sum2 += static_cast<i128>(sign) * s2 * t;
        s2 = t;

        std::swap(p, q);
        m = -m - 1;
        sign = -sign;
    }
}

u32 isqrt_u32(u32 n) {
    u32 r = static_cast<u32>(std::sqrt(static_cast<long double>(n)));
    while (static_cast<u64>(r + 1) * static_cast<u64>(r + 1) <= n) {
        ++r;
    }
    while (static_cast<u64>(r) * static_cast<u64>(r) > n) {
        --r;
    }
    return r;
}

i128 calc_smart(u32 val) {
    const u32 max_xy = isqrt_u32(val);
    i128 result1 = 0;
    i128 result2 = 0;

    for (u32 x = 2; x <= max_xy; ++x) {
        for (u32 y = 1; y < x; ++y) {
            const i64 x64 = x;
            const i64 y64 = y;
            const i64 v64 = val;

            const i64 split = v64 / (x64 + y64);
            result1 += static_cast<i128>(split) * (split - 1) / 2;

            const i64 s1 = v64 / x64 - split;
            const i64 s2 = (split <= max_xy) ? (static_cast<i64>(max_xy) - split) : 0;
            const FloorPair block = floor_sum_pair(-x64, y64, v64 - split * x64, s1, s2);
            result1 += block.sum1;

            if (split <= max_xy) {
                result2 += static_cast<i128>(split) * (split - 1) / 2;
                result2 += block.sum2;
            } else {
                result2 += static_cast<i128>(max_xy) * (max_xy - 1) / 2;
            }
        }
    }

    return 2 * result1 - result2;
}

void build_mobius(u32 limit, std::vector<i8>& mu, std::vector<i32>& prefix) {
    std::vector<u32> primes;
    std::vector<u32> lp(static_cast<std::size_t>(limit) + 1, 0);
    mu.assign(static_cast<std::size_t>(limit) + 1, 0);
    prefix.assign(static_cast<std::size_t>(limit) + 1, 0);

    mu[1] = 1;
    for (u32 i = 2; i <= limit; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes.push_back(i);
            mu[i] = -1;
        }
        for (const u32 p : primes) {
            const u64 v = static_cast<u64>(i) * p;
            if (v > limit) {
                break;
            }
            lp[static_cast<u32>(v)] = p;
            if (p == lp[i]) {
                mu[static_cast<u32>(v)] = 0;
                break;
            }
            mu[static_cast<u32>(v)] = static_cast<i8>(-mu[i]);
        }
    }

    for (u32 i = 1; i <= limit; ++i) {
        prefix[i] = prefix[i - 1] + static_cast<i32>(mu[i]);
    }
}

i128 base_term(u32 n) {
    const i64 m = n / 2;
    if ((n & 1U) == 0U) {
        return static_cast<i128>(3 * m - 2) * m;
    }
    return static_cast<i128>(3) * m * m + m;
}

u64 compute_s(u32 n,
              unsigned threads,
              u32 /*frontier_split_x*/,
              u32 /*chunk_divisor*/) {
    const u32 limit = n / 5;
    std::vector<i8> mu;
    std::vector<i32> prefix_mu;
    build_mobius(limit, mu, prefix_mu);

    std::vector<IntervalTask> tasks;
    tasks.reserve(static_cast<std::size_t>(2 * isqrt_u32(n) + 32));

    for (u32 l = 1; l <= limit;) {
        const u32 q = n / l;
        const u32 r = std::min<u32>(limit, n / q);
        const i32 mu_sum = prefix_mu[r] - prefix_mu[l - 1];
        if (mu_sum != 0 && q > 1) {
            tasks.push_back(IntervalTask{q, mu_sum});
        }
        l = r + 1;
    }

    i128 result = base_term(n);
    if (threads <= 1 || tasks.size() < 2) {
        for (const IntervalTask& task : tasks) {
            result += static_cast<i128>(2) * task.mu_sum * calc_smart(task.q);
        }
    } else {
        std::atomic<std::size_t> next_task(0);
        std::vector<i128> partial(static_cast<std::size_t>(threads), 0);
        std::vector<std::thread> pool;
        pool.reserve(static_cast<std::size_t>(threads));

        auto worker = [&](unsigned tid) {
            i128 local = 0;
            while (true) {
                const std::size_t idx = next_task.fetch_add(1, std::memory_order_relaxed);
                if (idx >= tasks.size()) {
                    break;
                }
                const IntervalTask& task = tasks[idx];
                local += static_cast<i128>(2) * task.mu_sum * calc_smart(task.q);
            }
            partial[static_cast<std::size_t>(tid)] = local;
        };

        for (unsigned t = 0; t < threads; ++t) {
            pool.emplace_back(worker, t);
        }
        for (auto& th : pool) {
            th.join();
        }
        for (const i128 v : partial) {
            result += v;
        }
    }

    const i128 total = static_cast<i128>(2) * result +
                       static_cast<i128>(n) * static_cast<i128>(n + 1) / 2;
    if (total < 0 || total > static_cast<i128>(std::numeric_limits<u64>::max())) {
        throw std::overflow_error("Result does not fit in u64");
    }
    return static_cast<u64>(total);
}

bool run_validation() {
    struct Check {
        u32 n;
        u64 expected;
    };

    const std::vector<Check> checks = {
        {1, 1},
        {10, 221},
        {100, 39'826},
        {1'000, 5'893'024},
    };

    bool ok = true;
    for (const Check& check : checks) {
        const u64 got = compute_s(check.n, 1, 128, 256);
        std::cout << "Validation S(" << check.n << ") = " << got;
        if (got == check.expected) {
            std::cout << " [PASS]";
        } else {
            std::cout << " [FAIL] expected " << check.expected;
            ok = false;
        }
        std::cout << '\n';
    }

    constexpr u32 kBruteN = 100;
    const u64 brute = brute_force_s(kBruteN);
    const u64 fast = compute_s(kBruteN, 1, 128, 256);
    std::cout << "Validation brute-vs-fast S(" << kBruteN << "): fast=" << fast
              << ", brute=" << brute;
    if (fast == brute) {
        std::cout << " [PASS]";
    } else {
        std::cout << " [FAIL]";
        ok = false;
    }
    std::cout << '\n';

    return ok;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_options(argc, argv, options)) {
        std::cerr << "Usage: ./Euler433 [--n=<positive-int>] [--threads=<positive-int>] "
                     "[--single-thread] [--no-validate] [--split=<positive-int>] "
                     "[--chunk-div=<positive-int>]\n";
        return 1;
    }

#if !defined(__OPTIMIZE__)
    std::cout << "Build warning: compiled without optimization flags. "
                 "Use -O3 for significantly faster runtime.\n";
#endif

    if (options.run_validation && !run_validation()) {
        return 1;
    }

    const unsigned threads = choose_thread_count(
        options.allow_multithreading, options.requested_threads, options.n / 2);

    std::cout << "Using " << threads << " thread(s)\n";
    if (threads > 1) {
        // Kept for CLI compatibility with previous scheduler-tuned versions.
        std::cout << "Scheduler: split_x=" << options.frontier_split_x
                  << ", chunk_div=" << options.chunk_divisor << '\n';
    }

    const auto start = std::chrono::steady_clock::now();
    const u64 result = compute_s(
        options.n, threads, options.frontier_split_x, options.chunk_divisor);
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = end - start;

    std::cout << "S(" << options.n << ") = " << result << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " seconds\n";
    return 0;
}
