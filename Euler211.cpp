#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Options {
    int limit = 64000000;
    bool run_checkpoints = true;
    int threads = 0;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 1 && options.threads >= 0;
}

inline bool is_perfect_square(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while (static_cast<unsigned __int128>(r + 1ULL) * static_cast<unsigned __int128>(r + 1ULL) <=
           x) {
        ++r;
    }
    while (static_cast<unsigned __int128>(r) * static_cast<unsigned __int128>(r) > x) {
        --r;
    }
    return static_cast<unsigned __int128>(r) * static_cast<unsigned __int128>(r) == x;
}

std::vector<u32> prime_list_upto(const int limit) {
    if (limit < 2) {
        return {};
    }
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;
    const int r = static_cast<int>(std::sqrt(static_cast<long double>(limit)));
    for (int p = 2; p <= r; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        const int step = p;
        int m = p * p;
        while (m <= limit) {
            is_prime[static_cast<std::size_t>(m)] = 0U;
            m += step;
        }
    }
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::max(1.0L, std::log((long double)limit))));
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] != 0U) {
            primes.push_back(static_cast<u32>(p));
        }
    }
    return primes;
}

u64 contribution(const std::vector<u32>& primes, const std::size_t from, const u64 s0,
                 const u64 n0, const u64 limit) {
    u64 result = 0ULL;
    const u64 max_n = limit - 1ULL;
    for (std::size_t i = from; i < primes.size(); ++i) {
        const u64 p = primes[i];
        if (n0 > max_n / p) {
            break;
        }

        const u64 p2 = p * p;
        u64 s = 1ULL;
        u64 pwr = p;

        while (n0 <= max_n / pwr) {
            const u64 n = n0 * pwr;
            s = s * p2 + 1ULL;
            const u64 s2 = s0 * s;

            if (is_perfect_square(s2)) {
                result += n;
            }
            result += contribution(primes, i + 1, s2, n, limit);

            if (pwr > max_n / p) {
                break;
            }
            pwr *= p;
        }
    }
    return result;
}

u64 contribution_root(const std::vector<u32>& primes, const std::size_t idx, const u64 limit) {
    const u64 max_n = limit - 1ULL;
    const u64 p = primes[idx];
    if (p > max_n) {
        return 0ULL;
    }

    const u64 p2 = p * p;
    u64 s = 1ULL;
    u64 pwr = p;
    u64 result = 0ULL;

    while (pwr <= max_n) {
        s = s * p2 + 1ULL;
        if (is_perfect_square(s)) {
            result += pwr;
        }
        result += contribution(primes, idx + 1, s, pwr, limit);
        if (pwr > max_n / p) {
            break;
        }
        pwr *= p;
    }

    return result;
}

struct WorkerArgs {
    const std::vector<u32>* primes = nullptr;
    u64 limit = 0ULL;
    std::atomic<std::size_t>* next_idx = nullptr;
    std::size_t max_idx = 0;
    std::size_t chunk = 1;
    u64 result = 0ULL;
};

void* worker_entry(void* raw) {
    auto* args = static_cast<WorkerArgs*>(raw);
    const auto& primes = *args->primes;
    u64 local = 0ULL;
    while (true) {
        const std::size_t start = args->next_idx->fetch_add(args->chunk, std::memory_order_relaxed);
        if (start >= args->max_idx) {
            break;
        }
        const std::size_t end = std::min(args->max_idx, start + args->chunk);
        for (std::size_t i = start; i < end; ++i) {
            local += contribution_root(primes, i, args->limit);
        }
    }
    args->result = local;
    return nullptr;
}

u64 solve(const int limit, int requested_threads) {
    if (limit <= 1) {
        return 0ULL;
    }

    const std::vector<u32> primes = prime_list_upto(limit - 1);
    u64 total = 1ULL;

    if (requested_threads <= 0) {
        long hw = sysconf(_SC_NPROCESSORS_ONLN);
        if (hw < 1) {
            hw = 1;
        }
        requested_threads = static_cast<int>(hw);
    }
    if (requested_threads < 1) {
        requested_threads = 1;
    }
    const std::size_t max_idx = primes.size();
    if (max_idx == 0) {
        return total;
    }

    std::size_t thread_count = static_cast<std::size_t>(requested_threads);
    if (thread_count > max_idx) {
        thread_count = max_idx;
    }
    if (thread_count <= 1) {
        for (std::size_t i = 0; i < max_idx; ++i) {
            total += contribution_root(primes, i, static_cast<u64>(limit));
        }
        return total;
    }

    std::atomic<std::size_t> next_idx(0);
    constexpr std::size_t chunk = 8;
    std::vector<pthread_t> tids(thread_count);
    std::vector<WorkerArgs> args(thread_count);

    bool create_failed = false;
    std::size_t started = 0;
    for (std::size_t t = 0; t < thread_count; ++t) {
        args[t].primes = &primes;
        args[t].limit = static_cast<u64>(limit);
        args[t].next_idx = &next_idx;
        args[t].max_idx = max_idx;
        args[t].chunk = chunk;
        args[t].result = 0ULL;

        if (pthread_create(&tids[t], nullptr, worker_entry, &args[t]) != 0) {
            create_failed = true;
            break;
        }
        ++started;
    }

    for (std::size_t t = 0; t < started; ++t) {
        pthread_join(tids[t], nullptr);
        total += args[t].result;
    }

    if (create_failed) {
        total = 1ULL;
        for (std::size_t i = 0; i < max_idx; ++i) {
            total += contribution_root(primes, i, static_cast<u64>(limit));
        }
    }

    return total;
}

bool run_checkpoints(const int threads) {
    if (solve(10, threads) != 1ULL) {
        std::cerr << "Checkpoint failed for limit=10" << '\n';
        return false;
    }
    if (solve(1000, threads) != 1304ULL) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
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
    if (options.run_checkpoints && !run_checkpoints(options.threads)) {
        return 2;
    }

    std::cout << solve(options.limit, options.threads) << '\n';
    return 0;
}
