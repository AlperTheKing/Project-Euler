#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

constexpr unsigned kDefaultN = 20U;
constexpr unsigned kMaxPracticalN = 24U;

struct Options {
    unsigned n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
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
        if (parse_unsigned_after_prefix(arg, "--n=", parsed_unsigned)) {
            options.n = parsed_unsigned;
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n < 4U) {
        std::cerr << "--n must be >= 4.\n";
        return false;
    }
    if (options.n > kMaxPracticalN) {
        std::cerr << "--n too large for this exact subset-enumeration implementation (max "
                  << kMaxPracticalN << ").\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const u64 workload_units) {
    constexpr u64 kMinUnitsForParallel = 200'000ULL;
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

long double ramvok_profit_real_cost(const int* values, const int k, const long double cost) {
    if (k <= 0) {
        return 0.0L;
    }

    const int max_value = values[k - 1];
    if (cost <= 0.0L) {
        return static_cast<long double>(max_value);
    }

    long double h_prev = 0.0L;
    long double best = 0.0L;

    for (int t = 1; t <= 10'000; ++t) {
        long double sum = 0.0L;
        for (int i = 0; i < k; ++i) {
            const long double x = static_cast<long double>(values[i]);
            sum += (x > h_prev) ? x : h_prev;
        }

        const long double h_cur = sum / static_cast<long double>(k);
        const long double candidate = h_cur - cost * static_cast<long double>(t);
        if (candidate > best) {
            best = candidate;
        }

        // Upper bound on any future profit at horizon >= t.
        const long double future_upper = static_cast<long double>(max_value) - cost * static_cast<long double>(t);
        if (future_upper <= best + 1e-16L) {
            break;
        }

        h_prev = h_cur;
    }

    return best;
}

void subset_best_profits_integer_costs(const std::uint32_t mask,
                                       const unsigned d,
                                       const unsigned max_cost,
                                       std::vector<long double>& output) {
    int values[24];
    int k = 0;

    for (unsigned bit = 0U; bit < d; ++bit) {
        if ((mask >> bit) & 1U) {
            values[k++] = static_cast<int>(bit + 1U);
        }
    }

    output.assign(max_cost + 1U, 0.0L);
    if (k == 0) {
        return;
    }

    const int max_value = values[k - 1];
    output[0] = static_cast<long double>(max_value);

    // For integer c>=1 and max roll <= d, t>d can never beat t=0 because
    // h_t <= max_value <= d, so h_t - c*t <= d - t < 0.
    const unsigned t_limit = d;

    long double h_prev = 0.0L;
    for (unsigned t = 1U; t <= t_limit; ++t) {
        long double sum = 0.0L;
        for (int i = 0; i < k; ++i) {
            const long double x = static_cast<long double>(values[i]);
            sum += (x > h_prev) ? x : h_prev;
        }

        const long double h_cur = sum / static_cast<long double>(k);
        for (unsigned c = 1U; c <= max_cost; ++c) {
            const long double candidate =
                h_cur - static_cast<long double>(c) * static_cast<long double>(t);
            if (candidate > output[c]) {
                output[c] = candidate;
            }
        }
        h_prev = h_cur;
    }
}

std::vector<long double> compute_average_reward_by_level(const unsigned d,
                                                         const unsigned max_cost,
                                                         const bool allow_multithreading,
                                                         const unsigned requested_threads) {
    const std::uint32_t total_masks = (1U << d);
    const u64 workload = static_cast<u64>(total_masks) - 1ULL;
    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, workload);

    const std::size_t stride = static_cast<std::size_t>(max_cost + 1U);
    const std::size_t table_size = static_cast<std::size_t>(d + 1U) * stride;

    std::vector<long double> global_sum(table_size, 0.0L);
    std::vector<u64> global_count(static_cast<std::size_t>(d + 1U), 0ULL);

    if (threads == 1U) {
        std::vector<long double> profits;
        for (std::uint32_t mask = 1U; mask < total_masks; ++mask) {
            const unsigned k = static_cast<unsigned>(__builtin_popcount(mask));
            subset_best_profits_integer_costs(mask, d, max_cost, profits);
            ++global_count[k];
            const std::size_t base = static_cast<std::size_t>(k) * stride;
            for (unsigned c = 0U; c <= max_cost; ++c) {
                global_sum[base + c] += profits[c];
            }
        }
    } else {
        std::vector<std::thread> workers;
        workers.reserve(threads);

        std::vector<std::vector<long double>> partial_sum(threads, std::vector<long double>(table_size, 0.0L));
        std::vector<std::vector<u64>> partial_count(threads, std::vector<u64>(d + 1U, 0ULL));

        for (unsigned t = 0U; t < threads; ++t) {
            const std::uint32_t begin =
                1U + static_cast<std::uint32_t>((static_cast<u64>(total_masks - 1U) * t) / threads);
            const std::uint32_t end =
                1U + static_cast<std::uint32_t>((static_cast<u64>(total_masks - 1U) * (t + 1U)) / threads);

            workers.emplace_back([&, t, begin, end]() {
                std::vector<long double> profits;
                for (std::uint32_t mask = begin; mask < end; ++mask) {
                    const unsigned k = static_cast<unsigned>(__builtin_popcount(mask));
                    subset_best_profits_integer_costs(mask, d, max_cost, profits);
                    ++partial_count[t][k];
                    const std::size_t base = static_cast<std::size_t>(k) * stride;
                    for (unsigned c = 0U; c <= max_cost; ++c) {
                        partial_sum[t][base + c] += profits[c];
                    }
                }
            });
        }

        for (auto& worker : workers) {
            worker.join();
        }

        for (unsigned t = 0U; t < threads; ++t) {
            for (unsigned k = 1U; k <= d; ++k) {
                global_count[k] += partial_count[t][k];
            }
            for (std::size_t i = 0; i < table_size; ++i) {
                global_sum[i] += partial_sum[t][i];
            }
        }
    }

    std::vector<long double> average(table_size, 0.0L);
    for (unsigned k = 1U; k <= d; ++k) {
        const u64 cnt = global_count[k];
        const long double inv_cnt = 1.0L / static_cast<long double>(cnt);
        const std::size_t base = static_cast<std::size_t>(k) * stride;
        for (unsigned c = 0U; c <= max_cost; ++c) {
            average[base + c] = global_sum[base + c] * inv_cnt;
        }
    }

    return average;
}

std::vector<long double> solve_super_by_cost(const unsigned d,
                                             const unsigned max_cost,
                                             const std::vector<long double>& average_reward_by_level) {
    const std::size_t stride = static_cast<std::size_t>(max_cost + 1U);
    std::vector<long double> s_by_cost(max_cost + 1U, 0.0L);

    std::vector<long double> lower(d, 0.0L);
    std::vector<long double> diag(d, 0.0L);
    std::vector<long double> upper(d, 0.0L);
    std::vector<long double> rhs(d, 0.0L);

    for (unsigned c = 0U; c <= max_cost; ++c) {
        for (unsigned k = 1U; k <= d; ++k) {
            const unsigned idx = k - 1U;
            const long double a_k =
                average_reward_by_level[static_cast<std::size_t>(k) * stride + c];

            if (k == d) {
                lower[idx] = -1.0L;
                diag[idx] = 1.0L;
                upper[idx] = 0.0L;
                rhs[idx] = a_k;
                continue;
            }

            const long double p = static_cast<long double>(k) / static_cast<long double>(d);
            const long double q = 1.0L - p;

            lower[idx] = (k == 1U) ? 0.0L : -p;
            diag[idx] = 1.0L;
            upper[idx] = -q;
            rhs[idx] = a_k;
        }

        for (unsigned i = 1U; i < d; ++i) {
            const long double factor = lower[i] / diag[i - 1U];
            diag[i] -= factor * upper[i - 1U];
            rhs[i] -= factor * rhs[i - 1U];
        }

        std::vector<long double> x(d, 0.0L);
        x[d - 1U] = rhs[d - 1U] / diag[d - 1U];
        for (int i = static_cast<int>(d) - 2; i >= 0; --i) {
            x[static_cast<std::size_t>(i)] =
                (rhs[static_cast<std::size_t>(i)] -
                 upper[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i + 1)]) /
                diag[static_cast<std::size_t>(i)];
        }

        s_by_cost[c] = x[d - 1U];
    }

    return s_by_cost;
}

bool nearly_equal(const long double a, const long double b, const long double eps) {
    return std::fabsl(a - b) <= eps;
}

bool run_checkpoints() {
    {
        const int values[] = {1, 2, 3, 4};
        const long double r = ramvok_profit_real_cost(values, 4, 0.2L);
        if (!nearly_equal(r, 2.65L, 1e-12L)) {
            std::cerr << "Checkpoint failed: R(4,0.2) expected 2.65, got "
                      << std::setprecision(18) << static_cast<double>(r) << '\n';
            return false;
        }
    }

    {
        constexpr unsigned d = 6U;
        constexpr unsigned max_cost = 1U;
        const std::vector<long double> avg =
            compute_average_reward_by_level(d, max_cost, false, 1U);
        const std::vector<long double> s = solve_super_by_cost(d, max_cost, avg);

        if (!nearly_equal(s[1], 208.3L, 1e-10L)) {
            std::cerr << "Checkpoint failed: S(6,1) expected 208.3, got "
                      << std::setprecision(18) << static_cast<double>(s[1]) << '\n';
            return false;
        }
    }

    return true;
}

long double compute_f(const unsigned n,
                      const bool allow_multithreading,
                      const unsigned requested_threads) {
    long double total = 0.0L;

    for (unsigned d = 4U; d <= n; ++d) {
        const std::vector<long double> avg =
            compute_average_reward_by_level(d, n, allow_multithreading, requested_threads);
        const std::vector<long double> s = solve_super_by_cost(d, n, avg);
        for (unsigned c = 0U; c <= n; ++c) {
            total += s[c];
        }
    }

    return total;
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

    const long double f_n = compute_f(options.n, options.allow_multithreading, options.requested_threads);
    const long long rounded = std::llround(f_n);

    std::cout << "F(" << options.n << ") = " << std::setprecision(18) << static_cast<double>(f_n)
              << '\n';
    std::cout << "Rounded: " << rounded << '\n';

    return 0;
}
