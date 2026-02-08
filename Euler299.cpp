#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 limit = 100000000ULL;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
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
            throw std::overflow_error("u64 parse overflow");
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        throw std::overflow_error("unsigned parse overflow");
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

unsigned pick_thread_count(const unsigned requested_threads) {
    if (requested_threads > 0U) {
        return requested_threads;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0U) {
        threads = 1U;
    }
    return threads;
}

template <typename Worker>
u64 parallel_sum_range(const int begin,
                       const int end,
                       unsigned threads,
                       const Worker& worker) {
    if (begin >= end) {
        return 0ULL;
    }

    const int work_items = end - begin;
    if (threads == 0U) {
        threads = 1U;
    }
    if (threads == 1U || work_items < 128) {
        return worker(begin, end);
    }

    threads = std::min<unsigned>(threads, static_cast<unsigned>(work_items));

    std::vector<std::thread> pool;
    pool.reserve(threads);

    std::vector<u64> partial(threads, 0ULL);

    for (unsigned tid = 0U; tid < threads; ++tid) {
        const int lo = begin + (work_items * static_cast<int>(tid)) /
                                   static_cast<int>(threads);
        const int hi = begin +
                       (work_items * static_cast<int>(tid + 1U)) /
                           static_cast<int>(threads);

        pool.emplace_back([&, tid, lo, hi]() { partial[tid] = worker(lo, hi); });
    }

    for (std::thread& th : pool) {
        th.join();
    }

    u64 total = 0ULL;
    for (const u64 x : partial) {
        total += x;
    }
    return total;
}

u64 count_family_x_eq_y(const u64 limit, unsigned threads) {
    if (limit <= 7ULL) {
        return 0ULL;
    }

    const u64 numerator = limit - 1ULL;
    const int m_max = static_cast<int>(std::sqrt(static_cast<long double>(limit))) + 3;

    auto worker = [limit, numerator](const int lo, const int hi) -> u64 {
        u64 local = 0ULL;

        for (int m = lo; m < hi; ++m) {
            const u64 m64 = static_cast<u64>(m);
            const u64 m2 = m64 * m64;

            // Orientation 1: u = g*m^2, v = 2*g*n^2 (m must be odd).
            if ((m & 1) != 0) {
                for (int n = 1;; ++n) {
                    const u64 n64 = static_cast<u64>(n);
                    const u64 s = m2 + 4ULL * m64 * n64 + 2ULL * n64 * n64;
                    if (s >= limit) {
                        break;
                    }
                    if (std::gcd(m, n) == 1) {
                        local += numerator / s;
                    }
                }
            }

            // Orientation 2: u = 2*g*m^2, v = g*n^2 (n must be odd).
            for (int n = 1;; ++n) {
                const u64 n64 = static_cast<u64>(n);
                const u64 s = 2ULL * m2 + 4ULL * m64 * n64 + n64 * n64;
                if (s >= limit) {
                    break;
                }
                if ((n & 1) != 0 && std::gcd(m, n) == 1) {
                    local += numerator / s;
                }
            }
        }

        return local;
    };

    return parallel_sum_range(1, m_max + 1, threads, worker);
}

u64 count_family_u_eq_v(const u64 limit, unsigned threads) {
    if (limit <= 13ULL) {
        return 0ULL;
    }

    const u64 numerator = limit - 1ULL;
    const int p_max = static_cast<int>(std::sqrt(static_cast<long double>(limit / 2ULL))) + 3;
    const int odd_count = (p_max + 1) / 2;

    auto worker = [limit, numerator](const int lo, const int hi) -> u64 {
        u64 local = 0ULL;

        for (int idx = lo; idx < hi; ++idx) {
            const int p = 2 * idx + 1;
            const u64 p64 = static_cast<u64>(p);
            const u64 p2 = p64 * p64;

            for (int q = 1;; ++q) {
                const u64 q64 = static_cast<u64>(q);
                const u64 s = p2 + 2ULL * p64 * q64 + 2ULL * q64 * q64;
                const u128 denom = static_cast<u128>(2ULL) * static_cast<u128>(s);
                if (denom >= static_cast<u128>(limit)) {
                    break;
                }
                if (std::gcd(p, q) == 1) {
                    local += numerator / static_cast<u64>(denom);
                }
            }
        }

        return local;
    };

    return parallel_sum_range(0, odd_count, threads, worker);
}

u64 count_triplets(const u64 limit, unsigned threads) {
    if (limit <= 0ULL) {
        return 0ULL;
    }

    if (limit < 1000000ULL) {
        threads = 1U;
    }

    const u64 family_x_eq_y = count_family_x_eq_y(limit, threads);
    const u64 family_u_eq_v = count_family_u_eq_v(limit, threads);
    return family_x_eq_y + family_u_eq_v;
}

void run_checkpoints(const unsigned threads) {
    struct Checkpoint {
        u64 limit;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {100ULL, 92ULL},
        {100000ULL, 320471ULL},
    };

    for (const Checkpoint& checkpoint : checkpoints) {
        const u64 got = count_triplets(checkpoint.limit, threads);
        if (got != checkpoint.expected) {
            throw std::runtime_error(
                "Checkpoint failed for limit=" + std::to_string(checkpoint.limit) +
                ": expected " + std::to_string(checkpoint.expected) +
                ", got " + std::to_string(got));
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }

        const unsigned threads = pick_thread_count(options.requested_threads);

        if (options.run_checkpoints) {
            run_checkpoints(threads);
        }

        const u64 answer = count_triplets(options.limit, threads);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
