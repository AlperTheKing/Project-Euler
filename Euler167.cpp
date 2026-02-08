#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Options {
    u64 target_k = 100000000000ULL;
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
        if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
            return false;
        }
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

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
            throw std::overflow_error("--k overflow");
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
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--k=", options.target_k)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

std::vector<u64> brute_ulam_terms(const u64 a, const u64 b, const std::size_t count) {
    if (count == 0U) {
        return {};
    }
    if (count == 1U) {
        return {a};
    }

    std::vector<u64> seq;
    seq.reserve(count);
    seq.push_back(a);
    seq.push_back(b);

    std::unordered_map<u64, u8> representation_count;
    representation_count.reserve(count * count / 2 + 64U);
    representation_count[a + b] = 1U;

    u64 candidate = b + 1U;

    while (seq.size() < count) {
        while (true) {
            const auto it = representation_count.find(candidate);
            const u8 ways = (it == representation_count.end()) ? 0U : it->second;
            if (ways == 1U) {
                break;
            }
            ++candidate;
        }

        const u64 next = candidate;
        for (const u64 x : seq) {
            if (x == next) {
                continue;
            }
            const u64 sum = x + next;
            u8& ways = representation_count[sum];
            if (ways < 2U) {
                ++ways;
            }
        }

        seq.push_back(next);
        candidate = next + 1U;
    }

    return seq;
}

std::vector<u64> brute_ulam_until_value(const u64 a, const u64 b, const u64 min_last_value) {
    std::vector<u64> seq;
    seq.reserve(256U);
    seq.push_back(a);
    seq.push_back(b);

    std::unordered_map<u64, u8> representation_count;
    representation_count.reserve(32768U);
    representation_count[a + b] = 1U;

    u64 candidate = b + 1U;

    while (seq.back() < min_last_value) {
        while (true) {
            const auto it = representation_count.find(candidate);
            const u8 ways = (it == representation_count.end()) ? 0U : it->second;
            if (ways == 1U) {
                break;
            }
            ++candidate;
        }

        const u64 next = candidate;
        for (const u64 x : seq) {
            if (x == next) {
                continue;
            }
            const u64 sum = x + next;
            u8& ways = representation_count[sum];
            if (ways < 2U) {
                ++ways;
            }
        }

        seq.push_back(next);
        candidate = next + 1U;
    }

    return seq;
}

struct UlamOddCycle {
    int v = 0;
    int gap = 0;
    u64 second_even = 0ULL;

    std::vector<u8> bits;
    std::vector<u64> prefix_ones;

    u64 cycle_t_start = 0ULL;
    u64 cycle_t_end = 0ULL;
    u64 cycle_len = 0ULL;
    u64 ones_before_cycle = 0ULL;
    u64 ones_per_cycle = 0ULL;

    std::vector<u64> cycle_prefix_ones;

    u64 odd_less_than_second_even = 0ULL;
    u64 second_even_index = 0ULL;

    explicit UlamOddCycle(const int vv) : v(vv), gap(vv + 1), second_even(2ULL * static_cast<u64>(vv + 1)) {
        if (v < 1) {
            throw std::invalid_argument("v must be positive");
        }
        if (gap <= 0 || gap > 30) {
            throw std::runtime_error("Unexpected recurrence gap");
        }

        const u64 seed_limit = 2ULL * static_cast<u64>(gap) + 1ULL;
        const std::vector<u64> seed = brute_ulam_until_value(2ULL, static_cast<u64>(v), seed_limit);

        std::unordered_set<u64> members;
        members.reserve(seed.size() * 2U + 16U);
        for (const u64 x : seed) {
            members.insert(x);
        }

        bits.resize(static_cast<std::size_t>(gap), 0U);
        for (int t = 0; t < gap; ++t) {
            const u64 odd_value = 2ULL * static_cast<u64>(t) + 1ULL;
            bits[static_cast<std::size_t>(t)] = (members.find(odd_value) != members.end()) ? 1U : 0U;
        }

        prefix_ones.reserve(1U << 14U);
        prefix_ones.push_back(0ULL);
        for (const u8 bit : bits) {
            prefix_ones.push_back(prefix_ones.back() + static_cast<u64>(bit));
        }

        odd_less_than_second_even = prefix_ones[static_cast<std::size_t>(gap)];
        second_even_index = 2ULL + odd_less_than_second_even;

        const u32 state_count = 1U << static_cast<u32>(gap);
        std::vector<int> seen(state_count, -1);

        u32 state = 0U;
        for (int i = 0; i < gap; ++i) {
            state = (state << 1U) | static_cast<u32>(bits[static_cast<std::size_t>(i)]);
        }

        seen[state] = 0;

        const u32 lower_mask = (gap == 1) ? 0U : ((1U << static_cast<u32>(gap - 1)) - 1U);

        int step = 0;
        int cycle_step_start = -1;
        int cycle_step_end = -1;

        while (true) {
            const u32 oldest = (state >> static_cast<u32>(gap - 1)) & 1U;
            const u32 newest = state & 1U;
            const u32 next_bit = oldest ^ newest;

            bits.push_back(static_cast<u8>(next_bit));
            prefix_ones.push_back(prefix_ones.back() + static_cast<u64>(next_bit));

            state = ((state & lower_mask) << 1U) | next_bit;
            ++step;

            const int prev = seen[state];
            if (prev >= 0) {
                cycle_step_start = prev;
                cycle_step_end = step;
                break;
            }
            seen[state] = step;
        }

        cycle_t_start = static_cast<u64>(gap - 1 + cycle_step_start);
        cycle_t_end = static_cast<u64>(gap - 1 + cycle_step_end);
        cycle_len = cycle_t_end - cycle_t_start;

        ones_before_cycle = prefix_ones[static_cast<std::size_t>(cycle_t_start + 1ULL)];
        ones_per_cycle = prefix_ones[static_cast<std::size_t>(cycle_t_end + 1ULL)] - ones_before_cycle;

        if (cycle_len == 0ULL || ones_per_cycle == 0ULL) {
            throw std::runtime_error("Degenerate cycle detected");
        }

        cycle_prefix_ones.assign(static_cast<std::size_t>(cycle_len + 1ULL), 0ULL);
        for (u64 i = 0ULL; i < cycle_len; ++i) {
            const u64 t = cycle_t_start + 1ULL + i;
            cycle_prefix_ones[static_cast<std::size_t>(i + 1ULL)] =
                cycle_prefix_ones[static_cast<std::size_t>(i)] + static_cast<u64>(bits[static_cast<std::size_t>(t)]);
        }
    }

    u64 odd_term_by_rank(const u64 odd_rank) const {
        if (odd_rank == 0ULL) {
            throw std::invalid_argument("odd rank must be positive");
        }

        u64 t = 0ULL;

        if (odd_rank <= ones_before_cycle) {
            const auto begin_it = prefix_ones.begin() + 1;
            const auto end_it = prefix_ones.begin() + static_cast<std::ptrdiff_t>(cycle_t_start + 2ULL);
            const auto it = std::lower_bound(begin_it, end_it, odd_rank);
            t = static_cast<u64>(std::distance(prefix_ones.begin(), it)) - 1ULL;
        } else {
            const u64 rem = odd_rank - ones_before_cycle;
            const u64 full_cycles = (rem - 1ULL) / ones_per_cycle;
            const u64 rem_inside_cycle = rem - full_cycles * ones_per_cycle;

            const auto begin_it = cycle_prefix_ones.begin() + 1;
            const auto end_it = cycle_prefix_ones.end();
            const auto it = std::lower_bound(begin_it, end_it, rem_inside_cycle);
            const u64 d = static_cast<u64>(std::distance(cycle_prefix_ones.begin(), it));

            t = cycle_t_start + full_cycles * cycle_len + d;
        }

        return 2ULL * t + 1ULL;
    }

    u64 kth_term(const u64 k) const {
        if (k == 0ULL) {
            throw std::invalid_argument("k must be positive");
        }
        if (k == 1ULL) {
            return 2ULL;
        }

        if (k == second_even_index) {
            return second_even;
        }
        if (k < second_even_index) {
            return odd_term_by_rank(k - 1ULL);
        }
        return odd_term_by_rank(k - 2ULL);
    }
};

bool run_checkpoints() {
    {
        const std::vector<u64> sample = brute_ulam_terms(1ULL, 2ULL, 7U);
        const std::vector<u64> expected{1ULL, 2ULL, 3ULL, 4ULL, 6ULL, 8ULL, 11ULL};
        if (sample != expected) {
            std::cerr << "Checkpoint failed: U(1,2) sample mismatch\n";
            return false;
        }
    }

    const std::vector<u64> ks{1ULL, 2ULL, 3ULL, 4ULL, 5ULL, 6ULL, 10ULL, 20ULL,
                              50ULL, 100ULL, 250ULL, 500ULL, 1000ULL, 1500ULL, 2000ULL};

    for (int n = 2; n <= 10; ++n) {
        const int v = 2 * n + 1;
        const UlamOddCycle fast(v);
        const std::vector<u64> brute = brute_ulam_terms(2ULL, static_cast<u64>(v), 2000U);

        for (const u64 k : ks) {
            const u64 a = fast.kth_term(k);
            const u64 b = brute[static_cast<std::size_t>(k - 1ULL)];
            if (a != b) {
                std::cerr << "Checkpoint failed: mismatch at v=" << v << ", k=" << k
                          << ", fast=" << a << ", brute=" << b << '\n';
                return false;
            }
        }

        std::unordered_set<u64> members;
        members.reserve(brute.size() * 2U + 16U);
        for (const u64 x : brute) {
            members.insert(x);
        }

        const u64 second_even = 2ULL * static_cast<u64>(v + 1);
        for (const u64 x : brute) {
            if ((x & 1ULL) == 0ULL && x != 2ULL && x != second_even) {
                std::cerr << "Checkpoint failed: unexpected extra even term for v=" << v
                          << ", value=" << x << '\n';
                return false;
            }
        }

        for (u64 odd = second_even + 1ULL; odd <= brute.back(); odd += 2ULL) {
            const bool lhs = (members.find(odd) != members.end());
            const bool rhs = (members.find(odd - 2ULL) != members.end()) ^
                             (members.find(odd - second_even) != members.end());
            if (lhs != rhs) {
                std::cerr << "Checkpoint failed: odd recurrence mismatch at v=" << v
                          << ", odd=" << odd << '\n';
                return false;
            }
        }
    }

    return true;
}

u64 solve_target_sum(const u64 target_k, const unsigned requested_threads) {
    constexpr int first_n = 2;
    constexpr int last_n = 10;
    constexpr int task_count = last_n - first_n + 1;

    std::vector<u64> values(task_count, 0ULL);

    const unsigned thread_count =
        std::max(1U, std::min<unsigned>(pick_thread_count(requested_threads), task_count));

    if (thread_count == 1U) {
        for (int i = 0; i < task_count; ++i) {
            const int n = first_n + i;
            const int v = 2 * n + 1;
            const UlamOddCycle solver(v);
            values[static_cast<std::size_t>(i)] = solver.kth_term(target_k);
        }
    } else {
        std::atomic<int> next_task(0);
        std::vector<std::thread> workers;
        workers.reserve(thread_count);

        for (unsigned t = 0; t < thread_count; ++t) {
            workers.emplace_back([&]() {
                while (true) {
                    const int task = next_task.fetch_add(1, std::memory_order_relaxed);
                    if (task >= task_count) {
                        return;
                    }

                    const int n = first_n + task;
                    const int v = 2 * n + 1;
                    const UlamOddCycle solver(v);
                    values[static_cast<std::size_t>(task)] = solver.kth_term(target_k);
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }
    }

    u64 sum = 0ULL;
    for (const u64 value : values) {
        sum += value;
    }
    return sum;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }

        if (options.run_checkpoints && !run_checkpoints()) {
            return 1;
        }

        const u64 answer = solve_target_sum(options.target_k, options.requested_threads);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
