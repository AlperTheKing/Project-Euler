#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetN = 1'053'779ULL;

struct Options {
    u64 n_limit = kTargetN;
    unsigned threads = 0;
    bool run_checkpoints = true;
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

    u64 parsed = 0;
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

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        u64 parsed = 0;
        if (parse_u64_after_prefix(arg, "--n=", parsed)) {
            options.n_limit = parsed;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--threads=", parsed)) {
            if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
                std::cerr << "Invalid thread count: " << parsed << "\n";
                return false;
            }
            options.threads = static_cast<unsigned>(parsed);
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        return false;
    }

    return true;
}

unsigned resolve_thread_count(unsigned requested) {
    if (requested > 0) {
        return requested;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }
    return threads;
}

inline bool within_radius_bound(const u64 multiplier,
                                const u64 q_value,
                                const u128 numerator_square) {
    const u128 scaled = static_cast<u128>(multiplier) * static_cast<u128>(q_value);
    return static_cast<u128>(3) * scaled * scaled <= numerator_square;
}

u64 floor_ratio_exact(const u64 numerator, const u64 q_value, const u128 numerator_square) {
    static const long double kSqrt3 = std::sqrt(3.0L);

    u64 estimate = static_cast<u64>(
        static_cast<long double>(numerator) /
        (static_cast<long double>(q_value) * kSqrt3));

    while (estimate < std::numeric_limits<u64>::max() &&
           within_radius_bound(estimate + 1, q_value, numerator_square)) {
        ++estimate;
    }
    while (!within_radius_bound(estimate, q_value, numerator_square)) {
        --estimate;
    }

    return estimate;
}

std::vector<u32> build_multiplier_table(const u64 numerator) {
    const u128 numerator_square = static_cast<u128>(numerator) * numerator;
    const u64 q_limit = floor_ratio_exact(numerator, 1, numerator_square);

    std::vector<u32> table(static_cast<std::size_t>(q_limit + 1), 0U);
    for (u64 q = 1; q <= q_limit; ++q) {
        table[static_cast<std::size_t>(q)] =
            static_cast<u32>(floor_ratio_exact(numerator, q, numerator_square));
    }
    return table;
}

template <typename Function>
u64 parallel_sum_over_range(const u64 first,
                            const u64 last,
                            unsigned threads,
                            const Function& contribution) {
    if (first > last) {
        return 0;
    }

    const u64 job_count = last - first + 1;
    if (threads == 0) {
        threads = 1;
    }
    threads = static_cast<unsigned>(std::min<u64>(threads, job_count));

    if (threads <= 1) {
        u64 total = 0;
        for (u64 i = first; i <= last; ++i) {
            total += contribution(i);
        }
        return total;
    }

    std::atomic<u64> next_index{first};
    std::vector<u64> partial(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            u64 local_sum = 0;
            while (true) {
                const u64 index = next_index.fetch_add(1, std::memory_order_relaxed);
                if (index > last) {
                    break;
                }
                local_sum += contribution(index);
            }
            partial[tid] = local_sum;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u64 total = 0;
    for (const u64 value : partial) {
        total += value;
    }
    return total;
}

u64 count_type_b(const std::vector<u32>& multiplier_table, const unsigned threads) {
    if (multiplier_table.size() <= 1) {
        return 0;
    }

    const u64 q_limit = static_cast<u64>(multiplier_table.size() - 1);

    u64 n_max = 0;
    for (u64 n = 1; n <= q_limit / (n + 1); ++n) {
        n_max = n;
    }

    auto contribution_for_n = [&](const u64 n) -> u64 {
        const u64 t_max = q_limit / n - n;
        u64 local = 0;

        for (u64 t = 1; t <= t_max; ++t) {
            if (t % 3 == 0) {
                continue;
            }
            if (std::gcd(n, t) != 1) {
                continue;
            }

            const u64 q = n * (n + t);  // q = m*n with m = n + t.
            local += multiplier_table[static_cast<std::size_t>(q)];
        }

        return local;
    };

    return parallel_sum_over_range(1, n_max, threads, contribution_for_n);
}

// Primitive 60-degree integer triangles come from primitive 120-degree triples:
// x = m^2 - n^2, y = 2mn + n^2, z = m^2 + mn + n^2 with gcd(m,n)=1 and (m-n)%3!=0.
// This yields two 60-degree families, parameterized with t = m-n:
// Type B: q = m*n = n*(n+t), r = k*sqrt(3)*q/2.
// Type A: q = (m-n)*(m+2n) = t*(3n+t), r = k*sqrt(3)*q/6.
// For fixed q, valid scales are k <= floor(C / (q*sqrt(3))) with C in {2N, 6N}.
u64 count_type_a(const std::vector<u32>& multiplier_table, const unsigned threads) {
    if (multiplier_table.size() <= 1) {
        return 0;
    }

    const u64 q_limit = static_cast<u64>(multiplier_table.size() - 1);

    u64 t_max = 0;
    for (u64 t = 1; t <= q_limit / (t + 3); ++t) {
        t_max = t;
    }

    auto contribution_for_t = [&](const u64 t) -> u64 {
        if (t % 3 == 0) {
            return 0;
        }

        const u64 n_max = (q_limit / t - t) / 3;
        u64 local = 0;

        for (u64 n = 1; n <= n_max; ++n) {
            if (std::gcd(n, t) != 1) {
                continue;
            }

            const u64 q = t * (3 * n + t);  // q = (m-n)*(m+2n).
            local += multiplier_table[static_cast<std::size_t>(q)];
        }

        return local;
    };

    return parallel_sum_over_range(1, t_max, threads, contribution_for_t);
}

u64 count_triangles(const u64 n_limit, const unsigned threads) {
    if (n_limit == 0) {
        return 0;
    }
    if (n_limit > std::numeric_limits<u64>::max() / 6ULL) {
        throw std::overflow_error("n is too large for safe arithmetic.");
    }

    const u64 numerator_b = 2ULL * n_limit;
    const std::vector<u32> table_b = build_multiplier_table(numerator_b);
    const u64 count_b = count_type_b(table_b, threads);

    const u64 numerator_a = 6ULL * n_limit;
    const std::vector<u32> table_a = build_multiplier_table(numerator_a);
    const u64 count_a = count_type_a(table_a, threads);

    return count_a + count_b;
}

bool run_validation_checkpoints() {
    struct Checkpoint {
        u64 n_limit;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {100, 1234},
        {1'000, 22'767},
        {10'000, 359'912},
    };

    for (const Checkpoint cp : checkpoints) {
        const u64 got = count_triangles(cp.n_limit, 1);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed for T(" << cp.n_limit << "): got " << got
                      << ", expected " << cp.expected << "\n";
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

    if (options.run_checkpoints && !run_validation_checkpoints()) {
        return 1;
    }

    const unsigned threads = resolve_thread_count(options.threads);
    const u64 answer = count_triangles(options.n_limit, threads);
    std::cout << answer << '\n';
    return 0;
}
