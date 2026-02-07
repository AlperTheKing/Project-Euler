#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

constexpr u64 kDefaultN = 1'000'000'000'000ULL;

struct Options {
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
    return true;
}

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while (static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) <= static_cast<u128>(x)) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > static_cast<u128>(x)) {
        --r;
    }
    return r;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }
    std::string out;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

// Computes sum_{m=L}^X floor(X/m).
u64 floor_sum_from(const u64 x, const u64 l) {
    if (l > x) {
        return 0ULL;
    }

    u128 total = 0U;
    if (static_cast<u128>(l) * static_cast<u128>(l) > static_cast<u128>(x)) {
        const u64 qmax = x / l;
        for (u64 q = 1ULL; q <= qmax; ++q) {
            u64 hi = x / q;
            u64 lo = x / (q + 1ULL) + 1ULL;
            if (hi < l) {
                continue;
            }
            if (lo < l) {
                lo = l;
            }
            if (lo <= hi) {
                total += static_cast<u128>(q) * static_cast<u128>(hi - lo + 1ULL);
            }
        }
    } else {
        for (u64 m = l; m <= x;) {
            const u64 q = x / m;
            const u64 r = x / q;
            total += static_cast<u128>(q) * static_cast<u128>(r - m + 1ULL);
            m = r + 1ULL;
        }
    }

    return static_cast<u64>(total);
}

struct MuDivCache {
    std::vector<std::vector<std::pair<u32, std::int8_t>>> divs;
};

MuDivCache build_mu_squarefree_divisors(const int nmax) {
    MuDivCache out;
    out.divs.resize(static_cast<std::size_t>(nmax) + 1ULL);
    out.divs[1].push_back({1U, 1});

    std::vector<int> spf(static_cast<std::size_t>(nmax) + 1ULL, 0);
    for (int i = 2; i <= nmax; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        spf[static_cast<std::size_t>(i)] = i;
        if (static_cast<u64>(i) * static_cast<u64>(i) > static_cast<u64>(nmax)) {
            continue;
        }
        for (u64 j = static_cast<u64>(i) * static_cast<u64>(i); j <= static_cast<u64>(nmax);
             j += static_cast<u64>(i)) {
            if (spf[static_cast<std::size_t>(j)] == 0) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    for (int i = 2; i <= nmax; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
        }
    }

    std::vector<int> primes;
    for (int x = 2; x <= nmax; ++x) {
        int y = x;
        primes.clear();
        while (y > 1) {
            const int p = spf[static_cast<std::size_t>(y)];
            primes.push_back(p);
            while (y % p == 0) {
                y /= p;
            }
        }

        const int k = static_cast<int>(primes.size());
        const int total_masks = 1 << k;
        auto& bucket = out.divs[static_cast<std::size_t>(x)];
        bucket.reserve(static_cast<std::size_t>(total_masks));

        for (int mask = 0; mask < total_masks; ++mask) {
            u32 d = 1U;
            int bit_count = 0;
            for (int i = 0; i < k; ++i) {
                if (((mask >> i) & 1) == 0) {
                    continue;
                }
                d *= static_cast<u32>(primes[static_cast<std::size_t>(i)]);
                ++bit_count;
            }
            const std::int8_t mu = static_cast<std::int8_t>((bit_count & 1) ? -1 : 1);
            bucket.push_back({d, mu});
        }
    }

    return out;
}

struct WorkerState {
    const MuDivCache* mu_cache = nullptr;
    u64 n = 0ULL;
    u64 odd_u_max = 0ULL;
    u64 even_w_max = 0ULL;
    std::atomic<u64> next_odd_u{1ULL};
    std::atomic<u64> next_even_w{1ULL};
};

i128 process_odd_u_block(const WorkerState& state, const u64 begin_u, const u64 end_u_exclusive) {
    i128 local = 0;
    for (u64 u = begin_u; u < end_u_exclusive; u += 2ULL) {
        if (u > state.odd_u_max) {
            break;
        }
        const auto& divisors = state.mu_cache->divs[static_cast<std::size_t>(u)];
        for (const auto& [d32, mu8] : divisors) {
            const u64 d = static_cast<u64>(d32);
            const u64 x = state.n / (u * d);
            if (x == 0ULL) {
                continue;
            }
            const u64 l = (u + d - 1ULL) / d; // ceil(u/d)
            const u64 contrib = floor_sum_from(x, l);
            if (mu8 > 0) {
                local += static_cast<i128>(contrib);
            } else {
                local -= static_cast<i128>(contrib);
            }
        }
    }
    return local;
}

i128 process_even_w_block(const WorkerState& state, const u64 begin_w, const u64 end_w_exclusive) {
    i128 local = 0;
    for (u64 w = begin_w; w < end_w_exclusive; ++w) {
        if (w > state.even_w_max) {
            break;
        }
        const auto& divisors = state.mu_cache->divs[static_cast<std::size_t>(w)];
        for (const auto& [d32, mu8] : divisors) {
            const u64 d = static_cast<u64>(d32);
            if ((d & 1ULL) == 0ULL) {
                continue;
            }

            const u64 x = state.n / (w * d);
            if (x == 0ULL) {
                continue;
            }
            const u64 l = (2ULL * w + d - 1ULL) / d; // ceil(2w/d)
            if (l > x) {
                continue;
            }

            // odd m only: sum_{odd m >= l} floor(x/m)
            const u64 odd_part =
                floor_sum_from(x, l) - floor_sum_from(x >> 1ULL, (l + 1ULL) >> 1ULL);
            if (mu8 > 0) {
                local += static_cast<i128>(odd_part);
            } else {
                local -= static_cast<i128>(odd_part);
            }
        }
    }
    return local;
}

u128 compute_f_sum(const u64 n, const MuDivCache& mu_cache, const unsigned threads) {
    const u64 odd_u_max = isqrt_u64(n);
    const u64 even_w_max = isqrt_u64(n / 2ULL);

    if (threads <= 1U) {
        const WorkerState state{&mu_cache, n, odd_u_max, even_w_max};
        i128 total = 0;
        total += process_odd_u_block(state, 1ULL, odd_u_max + 2ULL);
        total += process_even_w_block(state, 1ULL, even_w_max + 1ULL);
        return static_cast<u128>(total);
    }

    WorkerState state;
    state.mu_cache = &mu_cache;
    state.n = n;
    state.odd_u_max = odd_u_max;
    state.even_w_max = even_w_max;
    state.next_odd_u.store(1ULL, std::memory_order_relaxed);
    state.next_even_w.store(1ULL, std::memory_order_relaxed);

    constexpr u64 kOddBlockSize = 256ULL; // counts odd u's, stride 2
    constexpr u64 kEvenBlockSize = 512ULL;

    std::vector<std::thread> pool;
    std::vector<i128> partial(static_cast<std::size_t>(threads), 0);
    pool.reserve(static_cast<std::size_t>(threads));

    for (unsigned tid = 0U; tid < threads; ++tid) {
        pool.emplace_back([&, tid]() {
            i128 local = 0;
            while (true) {
                bool did_work = false;

                const u64 odd_start =
                    state.next_odd_u.fetch_add(2ULL * kOddBlockSize, std::memory_order_relaxed);
                if (odd_start <= state.odd_u_max) {
                    const u64 odd_end = odd_start + 2ULL * kOddBlockSize;
                    local += process_odd_u_block(state, odd_start, odd_end);
                    did_work = true;
                }

                const u64 even_start =
                    state.next_even_w.fetch_add(kEvenBlockSize, std::memory_order_relaxed);
                if (even_start <= state.even_w_max) {
                    const u64 even_end = even_start + kEvenBlockSize;
                    local += process_even_w_block(state, even_start, even_end);
                    did_work = true;
                }

                if (!did_work) {
                    break;
                }
            }
            partial[static_cast<std::size_t>(tid)] = local;
        });
    }

    for (std::thread& th : pool) {
        th.join();
    }

    i128 total = 0;
    for (const i128 value : partial) {
        total += value;
    }

    return static_cast<u128>(total);
}

u64 brute_f(const u64 n) {
    u64 total = 0ULL;
    for (u64 x = 1ULL; x <= n; ++x) {
        const u128 m = static_cast<u128>(2ULL) * static_cast<u128>(x) * static_cast<u128>(x);
        u64 c = 0ULL;
        for (u64 d = 1ULL; d <= x; ++d) {
            if (m % static_cast<u128>(d) == 0U) {
                ++c;
            }
        }
        total += c;
    }
    return total;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const u64 odd_u_max,
                             const u64 even_w_max) {
    if (!allow_multithreading) {
        return 1U;
    }

    const u64 workload = odd_u_max / 2ULL + even_w_max;
    if (workload < 200'000ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
        threads = std::min<unsigned>(threads, 16U);
    }
    return std::max(1U, threads);
}

bool run_checkpoints(const MuDivCache& mu_cache, const unsigned threads) {
    struct Checkpoint {
        u64 n;
        u64 expected;
    };
    const std::vector<Checkpoint> checkpoints = {
        {15ULL, 63ULL},
        {1'000ULL, 15'066ULL},
    };

    for (const Checkpoint cp : checkpoints) {
        const u128 got = compute_f_sum(cp.n, mu_cache, 1U);
        if (got != static_cast<u128>(cp.expected)) {
            std::cerr << "Checkpoint failed: F(" << cp.n << ") expected " << cp.expected
                      << ", got " << to_string_u128(got) << '\n';
            return false;
        }
    }

    const u64 brute_limit = 120ULL;
    const u64 brute_expected = brute_f(brute_limit);
    const u128 fast = compute_f_sum(brute_limit, mu_cache, 1U);
    if (fast != static_cast<u128>(brute_expected)) {
        std::cerr << "Brute-force cross-check failed at N=" << brute_limit
                  << ": expected " << brute_expected << ", got " << to_string_u128(fast) << '\n';
        return false;
    }

    if (threads > 1U) {
        const u64 thread_check_n = 10'000'000ULL;
        const u128 single = compute_f_sum(thread_check_n, mu_cache, 1U);
        const u128 multi = compute_f_sum(thread_check_n, mu_cache, threads);
        if (single != multi) {
            std::cerr << "Thread consistency failed at N=" << thread_check_n
                      << ": single=" << to_string_u128(single)
                      << ", multi=" << to_string_u128(multi) << '\n';
            return false;
        }
    }

    std::cout << "Validation checkpoints passed.\n";
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const u64 odd_u_max = isqrt_u64(options.n);
    const u64 even_w_max = isqrt_u64(options.n / 2ULL);
    const u64 limit_u64 = std::max(odd_u_max, even_w_max) + 8ULL;
    if (limit_u64 > static_cast<u64>(std::numeric_limits<int>::max())) {
        std::cerr << "N is too large for this implementation.\n";
        return 1;
    }
    const int divisor_cache_limit = static_cast<int>(limit_u64);
    const MuDivCache mu_cache = build_mu_squarefree_divisors(divisor_cache_limit);

    const unsigned threads = choose_thread_count(
        options.allow_multithreading, options.requested_threads, odd_u_max, even_w_max);

    if (options.run_checkpoints && !run_checkpoints(mu_cache, threads)) {
        return 1;
    }

    const u128 answer = compute_f_sum(options.n, mu_cache, threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
