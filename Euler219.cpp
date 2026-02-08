#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 n = 1'000'000'000ULL;
    u64 checkpoint_max_n = 200'000ULL;
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
                                 const std::string& prefix,
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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--checkpoint-max-n=", options.checkpoint_max_n)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n == 0ULL) {
        std::cerr << "--n must be >= 1.\n";
        return false;
    }

    return true;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }

    std::string digits;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

u128 choose_exact(unsigned n, unsigned k) {
    if (k > n) {
        return 0U;
    }
    k = std::min(k, n - k);

    u128 result = 1U;
    for (unsigned i = 1; i <= k; ++i) {
        result = (result * static_cast<u128>(n - k + i)) / static_cast<u128>(i);
    }

    return result;
}

std::vector<u128> build_cost_multiplicity(const unsigned max_cost) {
    std::vector<u128> count(static_cast<std::size_t>(max_cost) + 1ULL, 0U);
    count[0] = 1U;

    for (unsigned c = 1; c <= max_cost; ++c) {
        count[c] = count[c - 1];
        if (c >= 4U) {
            count[c] += count[c - 4U];
        }
    }

    return count;
}

u128 cost_multiplicity_combinatorial(const unsigned cost) {
    u128 total = 0U;

    for (unsigned ones = 0; 4U * ones <= cost; ++ones) {
        const unsigned zeros = cost - 4U * ones;
        const unsigned length = zeros + ones;
        total += choose_exact(length, ones);
    }

    return total;
}

u128 sum_smallest_internal_node_costs(const u64 internal_nodes) {
    if (internal_nodes == 0ULL) {
        return 0U;
    }

    std::vector<u128> multiplicity;
    multiplicity.reserve(128U);
    multiplicity.push_back(1U);

    u64 used = 0ULL;
    u128 sum = 0U;

    for (u64 cost = 0ULL; used < internal_nodes; ++cost) {
        if (cost >= multiplicity.size()) {
            u128 next = multiplicity[static_cast<std::size_t>(cost - 1ULL)];
            if (cost >= 4ULL) {
                next += multiplicity[static_cast<std::size_t>(cost - 4ULL)];
            }
            multiplicity.push_back(next);
        }

        const u64 remaining = internal_nodes - used;
        const u64 take = (multiplicity[static_cast<std::size_t>(cost)] >= static_cast<u128>(remaining))
                             ? remaining
                             : static_cast<u64>(multiplicity[static_cast<std::size_t>(cost)]);

        sum += static_cast<u128>(take) * static_cast<u128>(cost);
        used += take;
    }

    return sum;
}

u128 solve_cost(const u64 n) {
    if (n <= 1ULL) {
        return 0U;
    }

    const u64 internal_nodes = n - 1ULL;
    const u128 internal_sum = sum_smallest_internal_node_costs(internal_nodes);
    return static_cast<u128>(5ULL) * static_cast<u128>(internal_nodes) + internal_sum;
}

std::vector<u128> greedy_reference_costs(const u64 max_n) {
    std::vector<u128> costs(static_cast<std::size_t>(max_n) + 1ULL, 0U);
    if (max_n == 0ULL) {
        return costs;
    }

    std::priority_queue<u64, std::vector<u64>, std::greater<u64>> leaves;
    leaves.push(0ULL);

    costs[1] = 0U;
    u128 total = 0U;

    for (u64 n = 2ULL; n <= max_n; ++n) {
        const u64 c = leaves.top();
        leaves.pop();

        total += static_cast<u128>(c) + static_cast<u128>(5ULL);

        leaves.push(c + 1ULL);
        leaves.push(c + 4ULL);

        costs[static_cast<std::size_t>(n)] = total;
    }

    return costs;
}

unsigned pick_thread_count(const unsigned requested_threads, const u64 workload) {
    if (workload < 20'000ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 4U;
        }
    }

    if (threads == 0U) {
        threads = 1U;
    }

    if (threads > workload) {
        threads = static_cast<unsigned>(workload);
    }

    return std::max(1U, threads);
}

bool threaded_formula_vs_greedy_check(const std::vector<u128>& reference,
                                      const unsigned requested_threads) {
    if (reference.size() <= 1U) {
        return true;
    }

    const u64 max_n = static_cast<u64>(reference.size() - 1U);
    const unsigned threads = pick_thread_count(requested_threads, max_n);

    std::atomic<bool> ok(true);
    std::atomic<u64> mismatch_n(0ULL);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            for (u64 n = 1ULL + static_cast<u64>(t); n <= max_n; n += static_cast<u64>(threads)) {
                if (!ok.load(std::memory_order_relaxed)) {
                    return;
                }

                const u128 got = solve_cost(n);
                if (got != reference[static_cast<std::size_t>(n)]) {
                    ok.store(false, std::memory_order_relaxed);
                    mismatch_n.store(n, std::memory_order_relaxed);
                    return;
                }
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    if (!ok.load(std::memory_order_relaxed)) {
        const u64 n = mismatch_n.load(std::memory_order_relaxed);
        std::cerr << "Checkpoint failed at n=" << n
                  << ": formula=" << to_string_u128(solve_cost(n))
                  << ", greedy=" << to_string_u128(reference[static_cast<std::size_t>(n)])
                  << '\n';
        return false;
    }

    return true;
}

bool run_checkpoints(const Options& options) {
    if (solve_cost(6ULL) != 35U) {
        std::cerr << "Checkpoint failed: Cost(6) must be 35.\n";
        return false;
    }

    constexpr unsigned kCostRecurrenceCheckLimit = 70U;
    const std::vector<u128> recurrence = build_cost_multiplicity(kCostRecurrenceCheckLimit);
    for (unsigned c = 0; c <= kCostRecurrenceCheckLimit; ++c) {
        const u128 via_recurrence = recurrence[c];
        const u128 via_combinatorics = cost_multiplicity_combinatorial(c);
        if (via_recurrence != via_combinatorics) {
            std::cerr << "Checkpoint failed for cost=" << c
                      << ": recurrence=" << to_string_u128(via_recurrence)
                      << ", combinatorics=" << to_string_u128(via_combinatorics)
                      << '\n';
            return false;
        }
    }

    const u64 max_n = std::max<u64>(6ULL, options.checkpoint_max_n);
    const std::vector<u128> reference = greedy_reference_costs(max_n);

    if (reference[6] != 35U) {
        std::cerr << "Checkpoint failed: greedy Cost(6) mismatch.\n";
        return false;
    }

    if (!threaded_formula_vs_greedy_check(reference, options.requested_threads)) {
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

    try {
        if (options.run_checkpoints && !run_checkpoints(options)) {
            return 1;
        }

        const u128 answer = solve_cost(options.n);
        std::cout << "Cost(" << options.n << ") = " << to_string_u128(answer) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
