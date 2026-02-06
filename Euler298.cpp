#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <functional>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr int kDefaultAlphabet = 10;
constexpr int kDefaultCapacity = 5;
constexpr int kDefaultTurns = 50;
constexpr int kMaxCapacity = 5;
constexpr int kMaxLabelsEncodable = 15;

constexpr long double kMainExpected = 1.7688229423800932L;
constexpr long double kToleranceTight = 1e-15L;
constexpr long double kToleranceLoose = 1e-12L;
constexpr long double kMassTolerance = 5e-14L;

struct Options {
    int alphabet = kDefaultAlphabet;
    int capacity = kDefaultCapacity;
    int turns = kDefaultTurns;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct State {
    int l_size = 0;
    int r_size = 0;
    std::array<int, kMaxCapacity> l{};
    std::array<int, kMaxCapacity> r{};
};

struct CallResult {
    bool larry_hit = false;
    bool robin_hit = false;
    int delta = 0;  // Larry hit minus Robin hit.
};

struct Transition {
    int next_state = -1;
    int delta = 0;
    int weight = 0;
};

struct Graph {
    int alphabet = 0;
    int capacity = 0;
    std::vector<State> states;
    std::vector<std::vector<Transition>> transitions;
};

struct ExactSolveResult {
    long double expectation = 0.0L;
    long double probability_mass = 0.0L;
};

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
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
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
        if (parsed > static_cast<u64>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u32>(std::numeric_limits<int>::max())) {
        return false;
    }
    value = static_cast<int>(parsed);
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
        int parsed_int = 0;
        if (parse_int_after_prefix(arg, "--alphabet=", parsed_int)) {
            options.alphabet = parsed_int;
            continue;
        }
        if (parse_int_after_prefix(arg, "--capacity=", parsed_int)) {
            options.capacity = parsed_int;
            continue;
        }
        if (parse_int_after_prefix(arg, "--turns=", parsed_int)) {
            options.turns = parsed_int;
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

    if (options.alphabet < 1) {
        std::cerr << "--alphabet must be at least 1.\n";
        return false;
    }
    if (options.capacity < 1 || options.capacity > kMaxCapacity) {
        std::cerr << "--capacity must be in [1, " << kMaxCapacity << "].\n";
        return false;
    }
    if (options.capacity > options.alphabet) {
        std::cerr << "--capacity cannot exceed --alphabet.\n";
        return false;
    }
    if (options.turns < 0) {
        std::cerr << "--turns must be non-negative.\n";
        return false;
    }
    if (options.alphabet > kMaxLabelsEncodable) {
        std::cerr << "--alphabet must be <= " << kMaxLabelsEncodable
                  << " (encoding limit).\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload) {
    if (!allow_multithreading || workload < 50'000ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, 32U));
}

int find_index(const int* arr, const int size, const int value) {
    for (int i = 0; i < size; ++i) {
        if (arr[i] == value) {
            return i;
        }
    }
    return -1;
}

void canonicalize_state(State& state) {
    std::array<int, kMaxLabelsEncodable> remap{};
    remap.fill(-1);

    int next_id = 0;

    for (int i = 0; i < state.l_size; ++i) {
        const int x = state.l[static_cast<std::size_t>(i)];
        if (remap[static_cast<std::size_t>(x)] < 0) {
            remap[static_cast<std::size_t>(x)] = next_id++;
        }
    }
    for (int i = 0; i < state.r_size; ++i) {
        const int x = state.r[static_cast<std::size_t>(i)];
        if (remap[static_cast<std::size_t>(x)] < 0) {
            remap[static_cast<std::size_t>(x)] = next_id++;
        }
    }

    for (int i = 0; i < state.l_size; ++i) {
        const int x = state.l[static_cast<std::size_t>(i)];
        state.l[static_cast<std::size_t>(i)] = remap[static_cast<std::size_t>(x)];
    }
    for (int i = 0; i < state.r_size; ++i) {
        const int x = state.r[static_cast<std::size_t>(i)];
        state.r[static_cast<std::size_t>(i)] = remap[static_cast<std::size_t>(x)];
    }
}

u64 encode_state(const State& state, const int capacity) {
    u64 key = 0ULL;
    key |= static_cast<u64>(state.l_size);
    key |= (static_cast<u64>(state.r_size) << 3U);

    unsigned shift = 6U;
    for (int i = 0; i < capacity; ++i) {
        const u64 nibble = (i < state.l_size)
                               ? static_cast<u64>(state.l[static_cast<std::size_t>(i)] + 1)
                               : 0ULL;
        key |= (nibble << shift);
        shift += 4U;
    }

    for (int i = 0; i < capacity; ++i) {
        const u64 nibble = (i < state.r_size)
                               ? static_cast<u64>(state.r[static_cast<std::size_t>(i)] + 1)
                               : 0ULL;
        key |= (nibble << shift);
        shift += 4U;
    }

    return key;
}

bool same_state(const State& a, const State& b) {
    if (a.l_size != b.l_size || a.r_size != b.r_size) {
        return false;
    }

    for (int i = 0; i < a.l_size; ++i) {
        if (a.l[static_cast<std::size_t>(i)] != b.l[static_cast<std::size_t>(i)]) {
            return false;
        }
    }
    for (int i = 0; i < a.r_size; ++i) {
        if (a.r[static_cast<std::size_t>(i)] != b.r[static_cast<std::size_t>(i)]) {
            return false;
        }
    }
    return true;
}

bool has_contiguous_labels(const State& state) {
    std::array<unsigned char, kMaxLabelsEncodable> seen{};
    seen.fill(0U);

    int count = 0;
    int mx = -1;
    for (int i = 0; i < state.l_size; ++i) {
        const int x = state.l[static_cast<std::size_t>(i)];
        if (seen[static_cast<std::size_t>(x)] == 0U) {
            seen[static_cast<std::size_t>(x)] = 1U;
            ++count;
            mx = std::max(mx, x);
        }
    }
    for (int i = 0; i < state.r_size; ++i) {
        const int x = state.r[static_cast<std::size_t>(i)];
        if (seen[static_cast<std::size_t>(x)] == 0U) {
            seen[static_cast<std::size_t>(x)] = 1U;
            ++count;
            mx = std::max(mx, x);
        }
    }

    return (mx + 1) == count;
}

CallResult apply_call_in_place(State& state, const int called_label, const int capacity) {
    CallResult result;

    {
        const int idx = find_index(state.l.data(), state.l_size, called_label);
        if (idx >= 0) {
            result.larry_hit = true;
            ++result.delta;

            for (int i = idx; i + 1 < state.l_size; ++i) {
                state.l[static_cast<std::size_t>(i)] = state.l[static_cast<std::size_t>(i + 1)];
            }
            state.l[static_cast<std::size_t>(state.l_size - 1)] = called_label;
        } else {
            if (state.l_size < capacity) {
                state.l[static_cast<std::size_t>(state.l_size)] = called_label;
                ++state.l_size;
            } else {
                for (int i = 0; i + 1 < capacity; ++i) {
                    state.l[static_cast<std::size_t>(i)] = state.l[static_cast<std::size_t>(i + 1)];
                }
                state.l[static_cast<std::size_t>(capacity - 1)] = called_label;
            }
        }
    }

    {
        const int idx = find_index(state.r.data(), state.r_size, called_label);
        if (idx >= 0) {
            result.robin_hit = true;
            --result.delta;
        } else {
            if (state.r_size < capacity) {
                state.r[static_cast<std::size_t>(state.r_size)] = called_label;
                ++state.r_size;
            } else {
                for (int i = 0; i + 1 < capacity; ++i) {
                    state.r[static_cast<std::size_t>(i)] = state.r[static_cast<std::size_t>(i + 1)];
                }
                state.r[static_cast<std::size_t>(capacity - 1)] = called_label;
            }
        }
    }

    return result;
}

int intern_state(const State& state,
                 const int capacity,
                 std::unordered_map<u64, int>& id_by_key,
                 std::vector<State>& states,
                 std::vector<std::vector<Transition>>& transitions,
                 std::vector<unsigned char>& expanded,
                 bool* is_new = nullptr) {
    State canonical = state;
    canonicalize_state(canonical);
    const u64 key = encode_state(canonical, capacity);

    const auto it = id_by_key.find(key);
    if (it != id_by_key.end()) {
        const State& existing = states[static_cast<std::size_t>(it->second)];
        if (!same_state(existing, canonical)) {
            std::cerr << "Internal encoding collision detected.\n";
            std::abort();
        }
        if (is_new != nullptr) {
            *is_new = false;
        }
        return it->second;
    }

    const int id = static_cast<int>(states.size());
    id_by_key.emplace(key, id);
    states.push_back(canonical);
    transitions.emplace_back();
    expanded.push_back(0U);

    if (is_new != nullptr) {
        *is_new = true;
    }
    return id;
}

void expand_state(const int state_id,
                  const int alphabet,
                  const int capacity,
                  std::unordered_map<u64, int>& id_by_key,
                  std::vector<State>& states,
                  std::vector<std::vector<Transition>>& transitions,
                  std::vector<unsigned char>& expanded,
                  std::queue<int>& bfs_queue) {
    if (expanded[static_cast<std::size_t>(state_id)] != 0U) {
        return;
    }
    expanded[static_cast<std::size_t>(state_id)] = 1U;

    const State state = states[static_cast<std::size_t>(state_id)];

    std::array<unsigned char, kMaxLabelsEncodable> present{};
    present.fill(0U);

    for (int i = 0; i < state.l_size; ++i) {
        const int x = state.l[static_cast<std::size_t>(i)];
        present[static_cast<std::size_t>(x)] = 1U;
    }
    for (int i = 0; i < state.r_size; ++i) {
        const int x = state.r[static_cast<std::size_t>(i)];
        present[static_cast<std::size_t>(x)] = 1U;
    }

    std::vector<int> known_labels;
    known_labels.reserve(static_cast<std::size_t>(state.l_size + state.r_size));
    for (int x = 0; x < alphabet; ++x) {
        if (present[static_cast<std::size_t>(x)] != 0U) {
            known_labels.push_back(x);
        }
    }

    const int outside_count = alphabet - static_cast<int>(known_labels.size());

    std::vector<Transition> out;
    out.reserve(known_labels.size() + 1ULL);

    auto push_or_merge = [&out](const int next_id, const int delta, const int weight) {
        for (Transition& tr : out) {
            if (tr.next_state == next_id && tr.delta == delta) {
                tr.weight += weight;
                return;
            }
        }
        out.push_back(Transition{next_id, delta, weight});
    };

    for (const int called : known_labels) {
        State next = state;
        const CallResult result = apply_call_in_place(next, called, capacity);
        canonicalize_state(next);

        bool is_new = false;
        const int next_id =
            intern_state(next, capacity, id_by_key, states, transitions, expanded, &is_new);
        if (is_new) {
            bfs_queue.push(next_id);
        }

        push_or_merge(next_id, result.delta, 1);
    }

    if (outside_count > 0) {
        int representative = -1;
        for (int x = 0; x < alphabet; ++x) {
            if (present[static_cast<std::size_t>(x)] == 0U) {
                representative = x;
                break;
            }
        }

        State next = state;
        const CallResult result = apply_call_in_place(next, representative, capacity);
        canonicalize_state(next);

        bool is_new = false;
        const int next_id =
            intern_state(next, capacity, id_by_key, states, transitions, expanded, &is_new);
        if (is_new) {
            bfs_queue.push(next_id);
        }

        push_or_merge(next_id, result.delta, outside_count);
    }

    transitions[static_cast<std::size_t>(state_id)] = std::move(out);
}

Graph build_graph(const int alphabet, const int capacity) {
    Graph graph;
    graph.alphabet = alphabet;
    graph.capacity = capacity;

    std::unordered_map<u64, int> id_by_key;
    id_by_key.reserve(1024U);

    std::vector<unsigned char> expanded;
    expanded.reserve(1024U);

    State initial;
    (void)intern_state(initial,
                       capacity,
                       id_by_key,
                       graph.states,
                       graph.transitions,
                       expanded,
                       nullptr);

    std::queue<int> bfs_queue;
    bfs_queue.push(0);

    while (!bfs_queue.empty()) {
        const int id = bfs_queue.front();
        bfs_queue.pop();

        expand_state(id,
                     alphabet,
                     capacity,
                     id_by_key,
                     graph.states,
                     graph.transitions,
                     expanded,
                     bfs_queue);
    }

    return graph;
}

ExactSolveResult solve_exact(const Graph& graph,
                             const int turns,
                             const bool allow_multithreading,
                             const unsigned requested_threads) {
    const int state_count = static_cast<int>(graph.states.size());
    const int width = 2 * turns + 1;
    const int offset = turns;

    std::vector<long double> current(static_cast<std::size_t>(state_count) *
                                         static_cast<std::size_t>(width),
                                     0.0L);
    std::vector<long double> next(static_cast<std::size_t>(state_count) *
                                      static_cast<std::size_t>(width),
                                  0.0L);

    current[static_cast<std::size_t>(offset)] = 1.0L;

    for (int turn = 1; turn <= turns; ++turn) {
        std::fill(next.begin(), next.end(), 0.0L);

        const int d_prev_min = offset - (turn - 1);
        const int d_prev_max = offset + (turn - 1);

        const std::size_t workload = static_cast<std::size_t>(state_count) *
                                     static_cast<std::size_t>(2 * turn - 1) * 11ULL;
        const unsigned threads = choose_thread_count(allow_multithreading,
                                                     requested_threads,
                                                     workload);

        if (threads == 1U) {
            for (int sid = 0; sid < state_count; ++sid) {
                const std::size_t src_base =
                    static_cast<std::size_t>(sid) * static_cast<std::size_t>(width);

                for (const Transition& tr : graph.transitions[static_cast<std::size_t>(sid)]) {
                    const long double p = static_cast<long double>(tr.weight) /
                                          static_cast<long double>(graph.alphabet);
                    const std::size_t dst_base = static_cast<std::size_t>(tr.next_state) *
                                                 static_cast<std::size_t>(width);

                    for (int d = d_prev_min; d <= d_prev_max; ++d) {
                        const long double prob = current[src_base + static_cast<std::size_t>(d)];
                        if (prob == 0.0L) {
                            continue;
                        }

                        const int nd = d + tr.delta;
                        next[dst_base + static_cast<std::size_t>(nd)] += prob * p;
                    }
                }
            }
        } else {
            std::vector<std::vector<long double>> local_buffers(
                threads,
                std::vector<long double>(static_cast<std::size_t>(state_count) *
                                             static_cast<std::size_t>(width),
                                         0.0L));

            std::vector<std::thread> workers;
            workers.reserve(threads);

            for (unsigned t = 0U; t < threads; ++t) {
                const int begin = (state_count * static_cast<int>(t)) /
                                  static_cast<int>(threads);
                const int end = (state_count * static_cast<int>(t + 1U)) /
                                static_cast<int>(threads);

                workers.emplace_back([&, begin, end, t]() {
                    std::vector<long double>& local = local_buffers[t];

                    for (int sid = begin; sid < end; ++sid) {
                        const std::size_t src_base =
                            static_cast<std::size_t>(sid) * static_cast<std::size_t>(width);

                        for (const Transition& tr :
                             graph.transitions[static_cast<std::size_t>(sid)]) {
                            const long double p = static_cast<long double>(tr.weight) /
                                                  static_cast<long double>(graph.alphabet);
                            const std::size_t dst_base =
                                static_cast<std::size_t>(tr.next_state) *
                                static_cast<std::size_t>(width);

                            for (int d = d_prev_min; d <= d_prev_max; ++d) {
                                const long double prob =
                                    current[src_base + static_cast<std::size_t>(d)];
                                if (prob == 0.0L) {
                                    continue;
                                }

                                const int nd = d + tr.delta;
                                local[dst_base + static_cast<std::size_t>(nd)] += prob * p;
                            }
                        }
                    }
                });
            }

            for (std::thread& worker : workers) {
                worker.join();
            }

            for (std::size_t idx = 0ULL; idx < next.size(); ++idx) {
                long double acc = 0.0L;
                for (unsigned t = 0U; t < threads; ++t) {
                    acc += local_buffers[t][idx];
                }
                next[idx] = acc;
            }
        }

        current.swap(next);
    }

    long double expectation = 0.0L;
    long double mass = 0.0L;

    for (int sid = 0; sid < state_count; ++sid) {
        const std::size_t base = static_cast<std::size_t>(sid) * static_cast<std::size_t>(width);
        for (int i = 0; i < width; ++i) {
            const long double p = current[base + static_cast<std::size_t>(i)];
            if (p == 0.0L) {
                continue;
            }

            const int diff = i - offset;
            expectation += p * static_cast<long double>(std::abs(diff));
            mass += p;
        }
    }

    return ExactSolveResult{expectation, mass};
}

long double brute_force_expectation(const int alphabet, const int capacity, const int turns) {
    long double total_abs_diff = 0.0L;

    std::function<void(int, State&, int)> dfs =
        [&](const int step, State& state, const int diff) {
            if (step == turns) {
                total_abs_diff += static_cast<long double>(std::abs(diff));
                return;
            }

            for (int x = 0; x < alphabet; ++x) {
                State next = state;
                const CallResult result = apply_call_in_place(next, x, capacity);
                dfs(step + 1, next, diff + result.delta);
            }
        };

    State initial;
    dfs(0, initial, 0);

    long double total_sequences = 1.0L;
    for (int i = 0; i < turns; ++i) {
        total_sequences *= static_cast<long double>(alphabet);
    }

    return total_abs_diff / total_sequences;
}

bool compare_sets(const std::array<int, kMaxCapacity>& memory,
                  const int size,
                  std::initializer_list<int> expected_values) {
    std::vector<int> got;
    got.reserve(static_cast<std::size_t>(size));
    for (int i = 0; i < size; ++i) {
        got.push_back(memory[static_cast<std::size_t>(i)]);
    }
    std::sort(got.begin(), got.end());

    std::vector<int> expected(expected_values.begin(), expected_values.end());
    std::sort(expected.begin(), expected.end());

    return got == expected;
}

bool run_example_checkpoint() {
    State game;
    int cumulative_diff = 0;

    const std::vector<int> sequence = {1, 2, 4, 6, 1, 8, 10, 2, 4, 1};
    const std::vector<int> expected_larry = {0, 0, 0, 0, 1, 1, 1, 1, 1, 2};
    const std::vector<int> expected_robin = {0, 0, 0, 0, 1, 1, 1, 2, 3, 3};

    int larry_score = 0;
    int robin_score = 0;

    for (std::size_t i = 0ULL; i < sequence.size(); ++i) {
        const CallResult result = apply_call_in_place(game, sequence[i], kDefaultCapacity);
        cumulative_diff += result.delta;

        if (result.larry_hit) {
            ++larry_score;
        }
        if (result.robin_hit) {
            ++robin_score;
        }

        if (larry_score != expected_larry[i] || robin_score != expected_robin[i]) {
            std::cerr << "Example checkpoint mismatch on turn " << (i + 1ULL)
                      << ": got (L=" << larry_score << ", R=" << robin_score
                      << "), expected (L=" << expected_larry[i]
                      << ", R=" << expected_robin[i] << ").\n";
            return false;
        }
    }

    if (cumulative_diff != -1) {
        std::cerr << "Example checkpoint mismatch in total difference: got "
                  << cumulative_diff << ", expected -1.\n";
        return false;
    }

    if (!compare_sets(game.l, game.l_size, {1, 2, 4, 8, 10})) {
        std::cerr << "Example checkpoint mismatch in Larry's final memory set.\n";
        return false;
    }

    if (!compare_sets(game.r, game.r_size, {1, 4, 6, 8, 10})) {
        std::cerr << "Example checkpoint mismatch in Robin's final memory set.\n";
        return false;
    }

    return true;
}

bool run_bruteforce_cross_checkpoint() {
    const struct {
        int alphabet;
        int capacity;
        int turns;
    } cases[] = {
        {3, 2, 7},
        {4, 2, 8},
    };

    for (const auto& test_case : cases) {
        const Graph graph = build_graph(test_case.alphabet, test_case.capacity);
        const ExactSolveResult exact = solve_exact(graph, test_case.turns, false, 1U);
        const long double brute =
            brute_force_expectation(test_case.alphabet, test_case.capacity, test_case.turns);

        if (std::fabs(exact.expectation - brute) > 1e-15L) {
            std::cerr << "Brute-force cross-check failed for (alphabet="
                      << test_case.alphabet << ", capacity=" << test_case.capacity
                      << ", turns=" << test_case.turns << "): exact="
                      << std::setprecision(18) << exact.expectation
                      << ", brute=" << brute << '\n';
            return false;
        }

        if (std::fabs(exact.probability_mass - 1.0L) > kMassTolerance) {
            std::cerr << "Probability mass drift in brute-force cross-check case.\n";
            return false;
        }
    }

    return true;
}

bool run_thread_consistency_checkpoint(const Graph& graph) {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw < 2U) {
        return true;
    }

    const ExactSolveResult single = solve_exact(graph, kDefaultTurns, false, 1U);
    const ExactSolveResult threaded = solve_exact(graph, kDefaultTurns, true, 2U);

    if (std::fabs(single.expectation - threaded.expectation) > kToleranceLoose) {
        std::cerr << "Thread consistency mismatch: single=" << std::setprecision(18)
                  << single.expectation << ", threaded=" << threaded.expectation << '\n';
        return false;
    }

    return true;
}

bool run_graph_structure_checkpoint(const Graph& graph) {
    for (const State& state : graph.states) {
        if (!has_contiguous_labels(state)) {
            std::cerr << "Non-contiguous labels detected in canonical state.\n";
            return false;
        }
    }

    for (std::size_t sid = 0ULL; sid < graph.transitions.size(); ++sid) {
        int sum_w = 0;
        for (const Transition& tr : graph.transitions[sid]) {
            sum_w += tr.weight;
        }
        if (sum_w != graph.alphabet) {
            std::cerr << "Transition weight sum mismatch at state " << sid
                      << ": got " << sum_w << ", expected " << graph.alphabet << ".\n";
            return false;
        }
    }

    return true;
}

bool run_checkpoints(const Graph& graph) {
    if (!run_graph_structure_checkpoint(graph)) {
        std::cerr << "Graph-structure checkpoint failed.\n";
        return false;
    }

    if (!run_example_checkpoint()) {
        std::cerr << "Example checkpoint failed.\n";
        return false;
    }

    if (!run_bruteforce_cross_checkpoint()) {
        std::cerr << "Brute-force checkpoint failed.\n";
        return false;
    }

    if (!run_thread_consistency_checkpoint(graph)) {
        std::cerr << "Thread consistency checkpoint failed.\n";
        return false;
    }

    const ExactSolveResult main_exact = solve_exact(graph, kDefaultTurns, false, 1U);
    if (std::fabs(main_exact.expectation - kMainExpected) > kToleranceTight) {
        std::cerr << "Main-value checkpoint mismatch: got " << std::setprecision(18)
                  << main_exact.expectation << ", expected " << kMainExpected << '\n';
        return false;
    }

    if (std::fabs(main_exact.probability_mass - 1.0L) > kMassTolerance) {
        std::cerr << "Probability mass checkpoint failed for main case.\n";
        return false;
    }

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

    const Graph graph = build_graph(options.alphabet, options.capacity);

    if (options.run_checkpoints &&
        options.alphabet == kDefaultAlphabet &&
        options.capacity == kDefaultCapacity &&
        options.turns == kDefaultTurns) {
        if (!run_checkpoints(graph)) {
            return 1;
        }
    }

    const auto start = std::chrono::steady_clock::now();
    const ExactSolveResult result =
        solve_exact(graph, options.turns, options.allow_multithreading, options.requested_threads);
    const auto end = std::chrono::steady_clock::now();

    const std::chrono::duration<long double> elapsed = end - start;

    std::cout << std::fixed << std::setprecision(15)
              << "E(|L-R|) = " << result.expectation << '\n';
    std::cout << std::setprecision(8)
              << "Answer (8 d.p.): " << result.expectation << '\n';
    std::cout << std::setprecision(15)
              << "Probability mass: " << result.probability_mass << '\n';
    std::cout << "States in canonical graph: " << graph.states.size() << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " seconds\n";

    return 0;
}
