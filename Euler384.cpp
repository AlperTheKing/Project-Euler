#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr int kMinT = 2;
constexpr int kMaxT = 45;
constexpr u64 kCheckpointT = 54321ULL;
constexpr u64 kCheckpointC = 12345ULL;
constexpr u64 kCheckpointExpected = 1220847710ULL;
constexpr u64 kExpectedAnswer = 3354706415856332783ULL;

struct Options {
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct CountPair {
    u64 positive = 0ULL;
    u64 negative = 0ULL;
};

struct CacheKey {
    u64 t = 0ULL;
    i64 x = 0LL;

    bool operator==(const CacheKey& other) const {
        return t == other.t && x == other.x;
    }
};

struct CacheKeyHash {
    std::size_t operator()(const CacheKey& key) const {
        const std::size_t h1 = std::hash<u64>{}(key.t);
        const std::size_t h2 = std::hash<i64>{}(key.x);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
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

i64 floor_div(const i64 a, const i64 b) {
    if (a >= 0) {
        return a / b;
    }
    return -(((-a) + b - 1) / b);
}

int sign_for_b(const u64 n) {
    const unsigned parity = static_cast<unsigned>(__builtin_popcountll(n & (n >> 1ULL)) & 1U);
    return parity == 0U ? 1 : -1;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) {
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

class RudinShapiroSelector {
public:
    u64 g(const u64 t, const u64 c) {
        if (t == 0ULL || c == 0ULL || c > t) {
            return 0ULL;
        }

        i64 lo = 0;
        i64 hi = 1;

        while (count_total(t, hi) < c) {
            if (hi > std::numeric_limits<i64>::max() / 2LL) {
                hi = std::numeric_limits<i64>::max();
                break;
            }
            hi *= 2LL;
        }

        while (lo < hi) {
            const i64 mid = lo + (hi - lo) / 2LL;
            if (count_total(t, mid) >= c) {
                hi = mid;
            } else {
                lo = mid + 1LL;
            }
        }

        return static_cast<u64>(lo);
    }

private:
    CountPair prefix_counts(const u64 t, const i64 x) {
        if (t == 0ULL || x < 0LL) {
            return {0ULL, 0ULL};
        }

        if (t == 1ULL) {
            return {1ULL, 0ULL};
        }

        const CacheKey key{t, x};
        const auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        CountPair result;
        if ((t & 1ULL) == 0ULL) {
            const u64 u = t >> 1U;
            const i64 x1 = floor_div(x - 1LL, 4LL);
            const i64 x3 = floor_div(x - 3LL, 4LL);
            const CountPair a = prefix_counts(u, x1);
            const CountPair b = prefix_counts(u, x3);

            if ((u & 1ULL) == 0ULL) {
                result.positive = a.positive + b.positive;
                result.negative = a.negative + b.negative;
            } else {
                result.positive = a.positive + b.negative;
                result.negative = a.negative + b.positive;
            }
        } else {
            const u64 u = t >> 1U;
            const i64 x0 = floor_div(x, 4LL);
            const i64 x2 = floor_div(x - 2LL, 4LL);

            const CountPair u0 = prefix_counts(u, x0);
            const CountPair u2 = prefix_counts(u, x2);
            const CountPair v0 = prefix_counts(u + 1ULL, x0);
            const CountPair v2 = prefix_counts(u + 1ULL, x2);

            if ((u & 1ULL) == 1ULL) {
                // u odd.
                result.positive = u2.positive + v0.positive;
                result.negative = u0.negative + v2.positive;
            } else {
                // u even.
                result.positive = u2.negative + v0.positive;
                result.negative = u0.negative + v2.negative;
            }
        }

        memo_.emplace(key, result);
        return result;
    }

    u64 count_total(const u64 t, const i64 x) {
        const CountPair pair = prefix_counts(t, x);
        return pair.positive + pair.negative;
    }

    std::unordered_map<CacheKey, CountPair, CacheKeyHash> memo_;
};

std::vector<u64> fibonacci_values() {
    std::vector<u64> fib(static_cast<std::size_t>(kMaxT) + 1ULL, 0ULL);
    fib[0] = 1ULL;
    fib[1] = 1ULL;
    for (int i = 2; i <= kMaxT; ++i) {
        fib[static_cast<std::size_t>(i)] = fib[static_cast<std::size_t>(i - 1)] +
                                           fib[static_cast<std::size_t>(i - 2)];
    }
    return fib;
}

u64 solve_sum_gf_single_thread(const std::vector<u64>& fib) {
    RudinShapiroSelector selector;
    u64 sum = 0ULL;

    for (int t = kMinT; t <= kMaxT; ++t) {
        const u64 value = selector.g(fib[static_cast<std::size_t>(t)],
                                     fib[static_cast<std::size_t>(t - 1)]);
        sum += value;
    }

    return sum;
}

u64 solve_sum_gf_parallel(const std::vector<u64>& fib,
                          const bool allow_multithreading,
                          const unsigned requested_threads) {
    constexpr std::size_t kTaskCount = static_cast<std::size_t>(kMaxT - kMinT + 1);
    const unsigned thread_count = choose_thread_count(allow_multithreading,
                                                      requested_threads,
                                                      kTaskCount);
    if (thread_count == 1U) {
        return solve_sum_gf_single_thread(fib);
    }

    std::vector<u64> values(static_cast<std::size_t>(kMaxT) + 1ULL, 0ULL);
    std::atomic<int> next_t(kMinT);

    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (unsigned tid = 0U; tid < thread_count; ++tid) {
        workers.emplace_back([&]() {
            RudinShapiroSelector selector;
            while (true) {
                const int t = next_t.fetch_add(1);
                if (t > kMaxT) {
                    break;
                }

                const u64 value = selector.g(fib[static_cast<std::size_t>(t)],
                                             fib[static_cast<std::size_t>(t - 1)]);
                values[static_cast<std::size_t>(t)] = value;
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u64 sum = 0ULL;
    for (int t = kMinT; t <= kMaxT; ++t) {
        sum += values[static_cast<std::size_t>(t)];
    }
    return sum;
}

bool run_small_bruteforce_check() {
    constexpr int kMaxLevel = 18;
    std::vector<std::vector<u64>> positions(static_cast<std::size_t>(kMaxLevel) + 1ULL);
    std::vector<int> done(static_cast<std::size_t>(kMaxLevel) + 1ULL, 0);

    i64 current = 0;
    for (u64 n = 0ULL; n < 1'000'000ULL; ++n) {
        current += sign_for_b(n);
        if (current >= 1 && current <= kMaxLevel) {
            const std::size_t level = static_cast<std::size_t>(current);
            if (done[level] < static_cast<int>(level)) {
                positions[level].push_back(n);
                ++done[level];
            }
        }

        bool all_done = true;
        for (int level = 1; level <= kMaxLevel; ++level) {
            if (done[static_cast<std::size_t>(level)] < level) {
                all_done = false;
                break;
            }
        }
        if (all_done) {
            break;
        }
    }

    RudinShapiroSelector selector;
    for (int t = 1; t <= kMaxLevel; ++t) {
        for (int c = 1; c <= t; ++c) {
            const u64 got = selector.g(static_cast<u64>(t), static_cast<u64>(c));
            const u64 expected = positions[static_cast<std::size_t>(t)][static_cast<std::size_t>(c - 1)];
            if (got != expected) {
                std::cerr << "Small brute-force check failed at g(" << t << "," << c
                          << "): expected " << expected << ", got " << got << '\n';
                return false;
            }
        }
    }

    return true;
}

bool run_checkpoints(const Options& options) {
    RudinShapiroSelector selector;

    if (selector.g(3ULL, 3ULL) != 6ULL) {
        std::cerr << "Checkpoint failed: g(3,3) should be 6.\n";
        return false;
    }

    if (selector.g(4ULL, 2ULL) != 7ULL) {
        std::cerr << "Checkpoint failed: g(4,2) should be 7.\n";
        return false;
    }

    const u64 checkpoint = selector.g(kCheckpointT, kCheckpointC);
    if (checkpoint != kCheckpointExpected) {
        std::cerr << "Checkpoint failed: g(" << kCheckpointT << ',' << kCheckpointC << ") expected "
                  << kCheckpointExpected << ", got " << checkpoint << '\n';
        return false;
    }

    if (!run_small_bruteforce_check()) {
        return false;
    }

    const std::vector<u64> fib = fibonacci_values();
    const u64 single = solve_sum_gf_single_thread(fib);
    const u64 parallel = solve_sum_gf_parallel(fib,
                                               options.allow_multithreading,
                                               options.requested_threads);
    if (single != parallel) {
        std::cerr << "Parallel consistency check failed: single=" << single
                  << ", parallel=" << parallel << '\n';
        return false;
    }

    std::cout << "Checkpoints passed (threads="
              << choose_thread_count(options.allow_multithreading,
                                     options.requested_threads,
                                     static_cast<std::size_t>(kMaxT - kMinT + 1))
              << ").\n";

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

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const std::vector<u64> fib = fibonacci_values();
    const u64 answer = solve_sum_gf_parallel(fib,
                                             options.allow_multithreading,
                                             options.requested_threads);

    if (answer != kExpectedAnswer) {
        std::cerr << "Answer mismatch: expected " << kExpectedAnswer << ", got "
                  << answer << '\n';
        return 1;
    }

    std::cout << "Answer: " << answer << '\n';
    return 0;
}
