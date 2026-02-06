#include <algorithm>
#include <array>
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
using u128 = unsigned __int128;

constexpr int kFacesPerDie = 6;
constexpr int kThreshold = 19;  // Represents wins >= 19, i.e. probability > 1/2.
constexpr int kDefaultN = 30;

struct Options {
    int n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct Transition {
    std::uint8_t new_used_a = 0;
    std::uint8_t new_used_b = 0;
    std::uint8_t new_used_c = 0;
    std::uint8_t add_w_ba = 0;
    std::uint8_t add_w_cb = 0;
    std::uint8_t add_w_ac = 0;
};

using StateKey = std::uint32_t;
using StateMap = std::unordered_map<StateKey, u128>;

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

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            if (parsed_u64 > static_cast<u64>(std::numeric_limits<int>::max())) {
                std::cerr << "--n is too large for this implementation.\n";
                return false;
            }
            options.n = static_cast<int>(parsed_u64);
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units,
                             const std::size_t min_units_for_parallel) {
    if (!allow_multithreading || workload_units < min_units_for_parallel || workload_units < 2ULL) {
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

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value > 0) {
        const unsigned digit = static_cast<unsigned>(value % 10);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

constexpr int prefix_index(const int used_a, const int used_b, const int used_c) {
    return (used_a * 7 + used_b) * 7 + used_c;
}

constexpr StateKey pack_state(const int used_a,
                              const int used_b,
                              const int used_c,
                              const int w_ba,
                              const int w_cb,
                              const int w_ac) {
    return static_cast<StateKey>(used_a | (used_b << 3) | (used_c << 6) |
                                 (w_ba << 9) | (w_cb << 14) | (w_ac << 19));
}

inline void unpack_state(const StateKey key,
                         int& used_a,
                         int& used_b,
                         int& used_c,
                         int& w_ba,
                         int& w_cb,
                         int& w_ac) {
    used_a = key & 7;
    used_b = (key >> 3) & 7;
    used_c = (key >> 6) & 7;
    w_ba = (key >> 9) & 31;
    w_cb = (key >> 14) & 31;
    w_ac = (key >> 19) & 31;
}

std::vector<std::vector<Transition>> build_transitions() {
    std::vector<std::vector<Transition>> transitions(343);

    for (int used_a = 0; used_a <= kFacesPerDie; ++used_a) {
        for (int used_b = 0; used_b <= kFacesPerDie; ++used_b) {
            for (int used_c = 0; used_c <= kFacesPerDie; ++used_c) {
                auto& list = transitions[static_cast<std::size_t>(prefix_index(used_a, used_b, used_c))];

                const int rem_a = kFacesPerDie - used_a;
                const int rem_b = kFacesPerDie - used_b;
                const int rem_c = kFacesPerDie - used_c;

                for (int add_a = 0; add_a <= rem_a; ++add_a) {
                    const int new_used_a = used_a + add_a;
                    for (int add_b = 0; add_b <= rem_b; ++add_b) {
                        const int new_used_b = used_b + add_b;
                        for (int add_c = 0; add_c <= rem_c; ++add_c) {
                            const int new_used_c = used_c + add_c;

                            Transition tr;
                            tr.new_used_a = static_cast<std::uint8_t>(new_used_a);
                            tr.new_used_b = static_cast<std::uint8_t>(new_used_b);
                            tr.new_used_c = static_cast<std::uint8_t>(new_used_c);
                            tr.add_w_ba = static_cast<std::uint8_t>(add_b * used_a);
                            tr.add_w_cb = static_cast<std::uint8_t>(add_c * used_b);
                            tr.add_w_ac = static_cast<std::uint8_t>(add_a * used_c);
                            list.push_back(tr);
                        }
                    }
                }
            }
        }
    }

    return transitions;
}

inline void advance_single_state(const StateKey key,
                                 const u128 ways,
                                 const int values_left_after,
                                 const std::vector<std::vector<Transition>>& transitions,
                                 StateMap& out) {
    int used_a = 0;
    int used_b = 0;
    int used_c = 0;
    int w_ba = 0;
    int w_cb = 0;
    int w_ac = 0;
    unpack_state(key, used_a, used_b, used_c, w_ba, w_cb, w_ac);

    const auto& list = transitions[static_cast<std::size_t>(prefix_index(used_a, used_b, used_c))];
    for (const Transition& tr : list) {
        if (values_left_after == 0 &&
            (tr.new_used_a != kFacesPerDie || tr.new_used_b != kFacesPerDie || tr.new_used_c != kFacesPerDie)) {
            continue;
        }

        int new_w_ba = w_ba + tr.add_w_ba;
        int new_w_cb = w_cb + tr.add_w_cb;
        int new_w_ac = w_ac + tr.add_w_ac;

        if (new_w_ba > kThreshold) {
            new_w_ba = kThreshold;
        }
        if (new_w_cb > kThreshold) {
            new_w_cb = kThreshold;
        }
        if (new_w_ac > kThreshold) {
            new_w_ac = kThreshold;
        }

        // Remaining B/C/A faces can add at most 6 wins each for w_ba/w_cb/w_ac respectively.
        if (new_w_ba + (kFacesPerDie - tr.new_used_b) * kFacesPerDie < kThreshold) {
            continue;
        }
        if (new_w_cb + (kFacesPerDie - tr.new_used_c) * kFacesPerDie < kThreshold) {
            continue;
        }
        if (new_w_ac + (kFacesPerDie - tr.new_used_a) * kFacesPerDie < kThreshold) {
            continue;
        }

        const StateKey next_key =
            pack_state(tr.new_used_a, tr.new_used_b, tr.new_used_c, new_w_ba, new_w_cb, new_w_ac);
        out[next_key] += ways;
    }
}

void advance_layer(const StateMap& current,
                   const int values_left_after,
                   const std::vector<std::vector<Transition>>& transitions,
                   const bool allow_multithreading,
                   const unsigned requested_threads,
                   StateMap& next) {
    next.clear();

    constexpr std::size_t kMinUnitsForParallel = 35'000ULL;
    const unsigned thread_count = choose_thread_count(
        allow_multithreading, requested_threads, current.size(), kMinUnitsForParallel);

    if (thread_count == 1U) {
        next.reserve(current.size() * 8ULL + 256ULL);
        for (const auto& kv : current) {
            advance_single_state(kv.first, kv.second, values_left_after, transitions, next);
        }
        return;
    }

    std::vector<std::pair<StateKey, u128>> items;
    items.reserve(current.size());
    for (const auto& kv : current) {
        items.push_back(kv);
    }

    std::vector<StateMap> local_maps(thread_count);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    const std::size_t block_size = (items.size() + thread_count - 1ULL) / thread_count;
    for (unsigned tid = 0; tid < thread_count; ++tid) {
        workers.emplace_back([&, tid]() {
            const std::size_t begin = static_cast<std::size_t>(tid) * block_size;
            const std::size_t end = std::min(items.size(), begin + block_size);
            if (begin >= end) {
                return;
            }

            StateMap& local = local_maps[tid];
            local.reserve((end - begin) * 8ULL + 64ULL);

            for (std::size_t idx = begin; idx < end; ++idx) {
                advance_single_state(items[idx].first,
                                     items[idx].second,
                                     values_left_after,
                                     transitions,
                                     local);
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    std::size_t total_size = 0ULL;
    for (const StateMap& local : local_maps) {
        total_size += local.size();
    }
    next.reserve(total_size + 256ULL);

    for (const StateMap& local : local_maps) {
        for (const auto& kv : local) {
            next[kv.first] += kv.second;
        }
    }
}

u128 count_ordered_cycles_dp(const int n,
                             const bool allow_multithreading,
                             const unsigned requested_threads) {
    if (n <= 0) {
        return 0;
    }

    static const std::vector<std::vector<Transition>> transitions = build_transitions();

    StateMap current;
    StateMap next;
    current.reserve(4'096ULL);
    next.reserve(4'096ULL);
    current[pack_state(0, 0, 0, 0, 0, 0)] = 1;

    for (int value = 1; value <= n; ++value) {
        const int values_left_after = n - value;
        advance_layer(current,
                      values_left_after,
                      transitions,
                      allow_multithreading,
                      requested_threads,
                      next);
        current.swap(next);
    }

    const StateKey terminal = pack_state(kFacesPerDie, kFacesPerDie, kFacesPerDie,
                                         kThreshold, kThreshold, kThreshold);
    const auto it = current.find(terminal);
    return (it == current.end()) ? 0 : it->second;
}

u128 count_nontransitive_sets_dp(const int n,
                                 const bool allow_multithreading,
                                 const unsigned requested_threads) {
    const u128 ordered = count_ordered_cycles_dp(n, allow_multithreading, requested_threads);
    return ordered / 3;
}

void generate_dice_rec(const int n,
                       const int depth,
                       const int min_value,
                       std::array<std::uint8_t, kFacesPerDie>& current,
                       std::vector<std::array<std::uint8_t, kFacesPerDie>>& out) {
    if (depth == kFacesPerDie) {
        out.push_back(current);
        return;
    }

    for (int value = min_value; value <= n; ++value) {
        current[static_cast<std::size_t>(depth)] = static_cast<std::uint8_t>(value);
        generate_dice_rec(n, depth + 1, value, current, out);
    }
}

std::vector<std::array<std::uint8_t, kFacesPerDie>> generate_all_dice(const int n) {
    std::vector<std::array<std::uint8_t, kFacesPerDie>> dice;
    std::array<std::uint8_t, kFacesPerDie> current{};
    generate_dice_rec(n, 0, 1, current, dice);
    return dice;
}

inline bool bit_test(const std::vector<u64>& row, const int index) {
    return ((row[static_cast<std::size_t>(index >> 6)] >> (index & 63)) & 1ULL) != 0ULL;
}

u64 brute_force_nontransitive_sets(const int n,
                                   const bool allow_multithreading,
                                   const unsigned requested_threads) {
    if (n <= 0) {
        return 0ULL;
    }

    const std::vector<std::array<std::uint8_t, kFacesPerDie>> dice = generate_all_dice(n);
    const int m = static_cast<int>(dice.size());
    const int word_count = (m + 63) / 64;

    std::vector<std::vector<u64>> wins(static_cast<std::size_t>(m),
                                       std::vector<u64>(static_cast<std::size_t>(word_count), 0ULL));

    const unsigned pair_threads =
        choose_thread_count(allow_multithreading, requested_threads, static_cast<std::size_t>(m), 64ULL);

    {
        std::atomic<int> next_i(0);
        std::vector<std::thread> workers;
        workers.reserve(pair_threads);

        for (unsigned tid = 0; tid < pair_threads; ++tid) {
            workers.emplace_back([&]() {
                while (true) {
                    const int i = next_i.fetch_add(1, std::memory_order_relaxed);
                    if (i >= m) {
                        break;
                    }

                    auto& row = wins[static_cast<std::size_t>(i)];
                    for (int j = 0; j < m; ++j) {
                        if (i == j) {
                            continue;
                        }

                        int win_count = 0;
                        for (const std::uint8_t x : dice[static_cast<std::size_t>(i)]) {
                            for (const std::uint8_t y : dice[static_cast<std::size_t>(j)]) {
                                if (x > y) {
                                    ++win_count;
                                }
                            }
                        }

                        if (win_count >= kThreshold) {
                            row[static_cast<std::size_t>(j >> 6)] |= (1ULL << (j & 63));
                        }
                    }
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }
    }

    const unsigned triple_threads =
        choose_thread_count(allow_multithreading, requested_threads, static_cast<std::size_t>(m), 32ULL);
    std::vector<u64> partial(static_cast<std::size_t>(triple_threads), 0ULL);

    {
        std::atomic<int> next_i(0);
        std::vector<std::thread> workers;
        workers.reserve(triple_threads);

        for (unsigned tid = 0; tid < triple_threads; ++tid) {
            workers.emplace_back([&, tid]() {
                u64 local = 0ULL;

                while (true) {
                    const int i = next_i.fetch_add(1, std::memory_order_relaxed);
                    if (i >= m) {
                        break;
                    }

                    for (int j = i + 1; j < m; ++j) {
                        for (int k = j + 1; k < m; ++k) {
                            // For each die, one of the other two must beat it.
                            const bool i_beaten = bit_test(wins[static_cast<std::size_t>(j)], i) ||
                                                  bit_test(wins[static_cast<std::size_t>(k)], i);
                            if (!i_beaten) {
                                continue;
                            }

                            const bool j_beaten = bit_test(wins[static_cast<std::size_t>(i)], j) ||
                                                  bit_test(wins[static_cast<std::size_t>(k)], j);
                            if (!j_beaten) {
                                continue;
                            }

                            const bool k_beaten = bit_test(wins[static_cast<std::size_t>(i)], k) ||
                                                  bit_test(wins[static_cast<std::size_t>(j)], k);
                            if (k_beaten) {
                                ++local;
                            }
                        }
                    }
                }

                partial[static_cast<std::size_t>(tid)] = local;
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }
    }

    u64 total = 0ULL;
    for (const u64 count : partial) {
        total += count;
    }
    return total;
}

bool run_checkpoints(const Options& options) {
    const u128 sample = count_nontransitive_sets_dp(7,
                                                    options.allow_multithreading,
                                                    options.requested_threads);
    if (sample != static_cast<u128>(9'780ULL)) {
        std::cerr << "Checkpoint failed for N=7: expected 9780, got "
                  << to_string_u128(sample) << '\n';
        return false;
    }

    const u128 dp_n6 = count_nontransitive_sets_dp(6,
                                                   options.allow_multithreading,
                                                   options.requested_threads);
    const u64 brute_n6 = brute_force_nontransitive_sets(6,
                                                        options.allow_multithreading,
                                                        options.requested_threads);
    if (dp_n6 != static_cast<u128>(brute_n6)) {
        std::cerr << "Checkpoint failed for N=6: DP gives " << to_string_u128(dp_n6)
                  << ", brute force gives " << brute_n6 << '\n';
        return false;
    }

    const unsigned threads = choose_thread_count(options.allow_multithreading,
                                                 options.requested_threads,
                                                 std::size_t{500},
                                                 std::size_t{2});
    std::cout << "Checkpoints passed (threads=" << threads << ").\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 1;
    }

    const u128 answer = count_nontransitive_sets_dp(options.n,
                                                    options.allow_multithreading,
                                                    options.requested_threads);
    std::cout << "Answer: " << to_string_u128(answer) << '\n';
    return 0;
}
