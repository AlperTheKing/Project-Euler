#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <pthread.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
constexpr u32 MOD = 1'000'000'007;
constexpr unsigned MAX_SIDE = 11;
constexpr u32 NONE = std::numeric_limits<u32>::max();

void require(const bool condition, const char* description) {
    if (!condition) throw std::runtime_error(std::string("Check failed: ") + description);
}

u32 add(const u32 a, const u32 b) {
    const u32 sum = a + b;
    return sum >= MOD ? sum - MOD : sum;
}

u32 multiply(const u32 a, const u32 b) {
    return static_cast<u32>(u64{a} * b % MOD);
}

u32 power(u32 base, u32 exponent) {
    u32 result = 1;
    for (; exponent != 0; exponent >>= 1, base = multiply(base, base)) {
        if (exponent & 1) result = multiply(result, base);
    }
    return result;
}

u64 advance(const u64 state, const unsigned column, const unsigned side, const bool occupied) {
    std::array<unsigned, MAX_SIDE> labels{};
    for (unsigned i = 0; i < side; ++i) labels[i] = static_cast<unsigned>((state >> (4 * i)) & 15);
    if (!occupied) {
        labels[column] = 0;
    } else {
        const std::array<unsigned, 3> adjacent{
            column == 0 ? 0 : labels[column - 1], labels[column],
            column + 1 == side ? 0 : labels[column + 1]
        };
        unsigned merged = 15;
        for (const unsigned label : adjacent) {
            if (label != 0) merged = std::min(merged, label);
        }
        for (unsigned i = 0; i < side; ++i) {
            for (const unsigned label : adjacent) {
                if (label != 0 && labels[i] == label) labels[i] = merged;
            }
        }
        labels[column] = merged;
    }

    // Label 1 connects to the top edge; other component labels are canonicalized.
    std::array<unsigned, 16> canonical{0, 1};
    unsigned next = 2;
    bool reaches_top = false;
    u64 result = 0;
    for (unsigned i = 0; i < side; ++i) {
        reaches_top |= labels[i] == 1;
        if (labels[i] != 0 && canonical[labels[i]] == 0) canonical[labels[i]] = next++;
        result |= u64{canonical[labels[i]]} << (4 * i);
    }
    return reaches_top ? result : 0;
}

struct Incoming {
    std::vector<u32> offsets;
    std::vector<u32> sources;
};

struct Frontier {
    unsigned side;
    std::vector<u64> states;
    std::vector<std::pair<u32, unsigned>> initial;
    std::vector<std::array<u32, 2 * MAX_SIDE>> edges;
    std::vector<Incoming> incoming;

    explicit Frontier(const unsigned width) : side(width) {
        require(side >= 1 && side <= MAX_SIDE, "board side between 1 and 11");
        incoming.resize(side);
        std::unordered_map<u64, u32> ids;
        ids.reserve(std::size_t{1} << (side + 4));
        const auto intern = [&](const u64 state) {
            const auto entry = ids.emplace(state, static_cast<u32>(states.size()));
            if (entry.second) states.push_back(state);
            return entry.first->second;
        };
        for (unsigned mask = 1; mask < (1U << side); ++mask) {
            u64 state = 0;
            unsigned count = 0;
            for (unsigned column = 0; column < side; ++column) {
                if (mask & (1U << column)) {
                    state |= u64{1} << (4 * column);
                    ++count;
                }
            }
            initial.emplace_back(intern(state), count);
        }
        for (std::size_t i = 0; i < states.size(); ++i) {
            const u64 state = states[i];
            std::array<u32, 2 * MAX_SIDE> transitions{};
            for (unsigned column = 0; column < side; ++column) {
                for (unsigned occupied = 0; occupied < 2; ++occupied) {
                    const u64 next = advance(state, column, side, occupied != 0);
                    transitions[2 * column + occupied] = next == 0 ? NONE : intern(next);
                }
            }
            edges.push_back(transitions);
        }
        for (unsigned column = 0; column < side; ++column) {
            auto& reverse = incoming[column];
            reverse.offsets.resize(states.size() + 1);
            for (const auto& transitions : edges) {
                for (unsigned occupied = 0; occupied < 2; ++occupied) {
                    const u32 target = transitions[2 * column + occupied];
                    if (target != NONE) ++reverse.offsets[target + 1];
                }
            }
            for (std::size_t i = 0; i < states.size(); ++i) reverse.offsets[i + 1] += reverse.offsets[i];
            reverse.sources.resize(reverse.offsets.back());
            auto positions = reverse.offsets;
            for (u32 source = 0; source < states.size(); ++source) {
                for (unsigned occupied = 0; occupied < 2; ++occupied) {
                    const u32 target = edges[source][2 * column + occupied];
                    if (target != NONE) reverse.sources[positions[target]++] = source;
                }
            }
        }
    }
};

struct Barrier {
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t changed = PTHREAD_COND_INITIALIZER;
    unsigned arrived = 0;
    unsigned generation = 0;

    ~Barrier() {
        pthread_cond_destroy(&changed);
        pthread_mutex_destroy(&mutex);
    }

    void wait(const unsigned participants) {
        if (participants == 1) return;
        pthread_mutex_lock(&mutex);
        const unsigned previous = generation;
        if (++arrived == participants) {
            arrived = 0;
            ++generation;
            pthread_cond_broadcast(&changed);
        } else {
            while (generation == previous) pthread_cond_wait(&changed, &mutex);
        }
        pthread_mutex_unlock(&mutex);
    }
};

struct Work {
    const Frontier& frontier;
    u32* values;
    u32* scratch;
    unsigned stride;
    unsigned threads;
    Barrier barrier;

    Work(const Frontier& graph, u32* first, u32* second, const unsigned length, const unsigned count)
        : frontier(graph), values(first), scratch(second), stride(length), threads(count) {}
};

struct Task {
    Work* work;
    unsigned index;
};

void* count_worker(void* argument) {
    const auto& task = *static_cast<Task*>(argument);
    auto& work = *task.work;
    pthread_mutex_lock(&work.barrier.mutex);
    pthread_mutex_unlock(&work.barrier.mutex);
    const auto& graph = work.frontier;
    const std::size_t first = graph.states.size() * task.index / work.threads;
    const std::size_t last = graph.states.size() * (task.index + 1) / work.threads;
    u32* current = work.values;
    u32* next = work.scratch;
    for (unsigned cell = graph.side; cell < graph.side * graph.side; ++cell) {
        const unsigned column = cell % graph.side;
        const unsigned upper = std::min(work.stride - 1, cell + 1);
        const auto& incoming = graph.incoming[column];
        for (std::size_t target = first; target < last; ++target) {
            const unsigned occupied = ((graph.states[target] >> (4 * column)) & 15) != 0;
            u32* output = next + target * work.stride;
            std::fill(output, output + upper + 1, 0);
            for (u32 j = incoming.offsets[target]; j < incoming.offsets[target + 1]; ++j) {
                const u32* input = current + std::size_t{incoming.sources[j]} * work.stride;
                for (unsigned degree = occupied; degree <= upper; ++degree) {
                    output[degree] = add(output[degree], input[degree - occupied]);
                }
            }
        }
        work.barrier.wait(work.threads);
        std::swap(current, next);
    }
    return nullptr;
}

std::vector<u32> crossing_counts(const Frontier& graph, unsigned threads) {
    const unsigned stride = graph.side * graph.side / 2 + 1;
    std::vector<u32> values(graph.states.size() * stride), scratch(values.size());
    for (const auto& entry : graph.initial) {
        if (entry.second < stride) values[std::size_t{entry.first} * stride + entry.second] = 1;
    }
    if (threads == 0) {
        threads = graph.side < 9 ? 1 : std::min(8U, std::max(1U, std::thread::hardware_concurrency()));
    }
    threads = static_cast<unsigned>(std::min<std::size_t>(threads, graph.states.size()));
    Work work(graph, values.data(), scratch.data(), stride, threads);
    std::vector<Task> tasks;
    for (unsigned i = 0; i < threads; ++i) tasks.push_back({&work, i});
    std::vector<pthread_t> workers(threads - 1);
    pthread_mutex_lock(&work.barrier.mutex);
    unsigned started = 1;
    for (; started < threads; ++started) {
        if (pthread_create(&workers[started - 1], nullptr, count_worker, &tasks[started]) != 0) break;
    }
    work.threads = started;
    pthread_mutex_unlock(&work.barrier.mutex);
    count_worker(&tasks[0]);
    bool joined = true;
    for (unsigned i = 1; i < started; ++i) {
        if (pthread_join(workers[i - 1], nullptr) != 0) joined = false;
    }
    require(joined, "pthread join");
    std::vector<u32> result(stride);
    for (std::size_t state = 0; state < graph.states.size(); ++state) {
        for (unsigned degree = 0; degree < stride; ++degree) {
            result[degree] = add(result[degree], values[state * stride + degree]);
        }
    }
    return result;
}

u32 game_count(const unsigned cells, const std::vector<u32>& crossings) {
    std::vector<u32> factorial(cells + 1, 1), inverse(cells + 1, 1), probability(cells / 2 + 2);
    for (unsigned i = 1; i <= cells; ++i) factorial[i] = multiply(factorial[i - 1], i);
    inverse[cells] = power(factorial[cells], MOD - 2);
    for (unsigned i = cells; i != 0; --i) inverse[i - 1] = multiply(inverse[i], i);
    for (unsigned k = 0; k <= cells / 2; ++k) {
        probability[k] = multiply(crossings[k], multiply(inverse[cells],
            multiply(factorial[k], factorial[cells - k])));
    }
    if (cells % 2 != 0) probability[cells / 2 + 1] = (1 + MOD - probability[cells / 2]) % MOD;

    // With p_k the crossing probability, wins on moves 2k+1 and 2k+2 total
    // (cells)_(2k+1) * (cells - 2k) * (p_(k+1) - p_k).
    u32 result = 0;
    for (unsigned k = 0; k < (cells + 1) / 2; ++k) {
        const u32 difference = (probability[k + 1] + MOD - probability[k]) % MOD;
        const u32 weight = multiply(multiply(factorial[cells], inverse[cells - 2 * k - 1]), cells - 2 * k);
        result = add(result, multiply(weight, difference));
    }
    return result;
}

struct SmallBoard {
    unsigned side;
    unsigned cells;
    std::vector<u32> adjacent;
    std::array<u32, 2> first{}, last{};

    explicit SmallBoard(const unsigned n) : side(n), cells(n * n), adjacent(cells) {
        constexpr std::array<int, 6> DR{-1, -1, 0, 0, 1, 1};
        constexpr std::array<int, 6> DC{0, 1, -1, 1, -1, 0};
        for (unsigned cell = 0; cell < cells; ++cell) {
            const unsigned row = cell / side, column = cell % side;
            const u32 bit = u32{1} << cell;
            if (row == 0) first[0] |= bit;
            if (row + 1 == side) last[0] |= bit;
            if (column == 0) first[1] |= bit;
            if (column + 1 == side) last[1] |= bit;
            for (unsigned d = 0; d < DR.size(); ++d) {
                const int r = static_cast<int>(row) + DR[d], c = static_cast<int>(column) + DC[d];
                if (r >= 0 && c >= 0 && r < static_cast<int>(side) && c < static_cast<int>(side)) {
                    adjacent[cell] |= u32{1} << (static_cast<unsigned>(r) * side + static_cast<unsigned>(c));
                }
            }
        }
    }

    bool wins(const u32 stones, const unsigned player) const {
        u32 connected = stones & first[player];
        for (;;) {
            if (connected & last[player]) return true;
            u32 next = connected;
            for (unsigned cell = 0; cell < cells; ++cell) {
                if (connected & (u32{1} << cell)) next |= adjacent[cell] & stones;
            }
            if (next == connected) return false;
            connected = next;
        }
    }

    u64 games(const u32 white, const u32 black, const unsigned player, std::vector<u64>& memo) const {
        u64& result = memo[white | (black << cells)];
        if (result != std::numeric_limits<u64>::max()) return result;
        if (wins(white, 0) || wins(black, 1)) return result = 1;
        result = 0;
        for (unsigned cell = 0; cell < cells; ++cell) {
            const u32 bit = u32{1} << cell;
            if ((white | black) & bit) continue;
            result += player == 0 ? games(white | bit, black, 1, memo) : games(white, black | bit, 0, memo);
        }
        return result;
    }
};

u32 evaluate_crossings(const Frontier& graph, const u32 weight) {
    std::vector<u32> current(graph.states.size()), next(current.size());
    for (const auto& entry : graph.initial) current[entry.first] = power(weight, entry.second);
    for (unsigned cell = graph.side; cell < graph.side * graph.side; ++cell) {
        std::fill(next.begin(), next.end(), 0);
        const unsigned column = cell % graph.side;
        for (std::size_t source = 0; source < current.size(); ++source) {
            for (unsigned occupied = 0; occupied < 2; ++occupied) {
                const u32 target = graph.edges[source][2 * column + occupied];
                if (target == NONE) continue;
                next[target] = add(next[target], occupied == 0 ? current[source] : multiply(current[source], weight));
            }
        }
        current.swap(next);
    }
    u32 result = 0;
    for (const u32 value : current) result = add(result, value);
    return result;
}

void run_tests() {
    constexpr std::array<u64, 4> GIVEN{1, 24, 257'760, 10'663'950'528'000};
    for (unsigned side = 1; side <= 4; ++side) {
        const Frontier graph(side);
        const auto actual = crossing_counts(graph, 1);
        const SmallBoard board(side);
        std::vector<u32> expected(board.cells + 1);
        const u32 all = (u32{1} << board.cells) - 1;
        for (u32 stones = 0; stones <= all; ++stones) {
            const bool white_wins = board.wins(stones, 0);
            require(white_wins != board.wins(all ^ stones, 1), "exactly one winner on a full board");
            unsigned count = 0;
            for (u32 remaining = stones; remaining != 0; remaining &= remaining - 1) ++count;
            if (white_wins) ++expected[count];
        }
        require(std::equal(actual.begin(), actual.end(), expected.begin()), "all subsets through side 4");
        const u32 answer = game_count(board.cells, actual);
        require(answer == GIVEN[side - 1] % MOD, "given game count");
        if (side <= 3) {
            std::vector<u64> memo(std::size_t{1} << (2 * board.cells), std::numeric_limits<u64>::max());
            require(answer == board.games(0, 0, 0, memo) % MOD, "direct game tree through side 3");
        }
    }

    const Frontier graph(MAX_SIDE);
    const auto serial = crossing_counts(graph, 1);
    require(serial == crossing_counts(graph, 4), "serial and pthread agreement");
    for (unsigned k = 0; k < MAX_SIDE; ++k) require(serial[k] == 0, "at least one stone in every row");
    require(serial[MAX_SIDE] == (MAX_SIDE + 1) * (1U << (MAX_SIDE - 2)), "minimal crossing paths");

    const unsigned cells = MAX_SIDE * MAX_SIDE;
    std::vector<u32> binomial(cells + 1), full(cells + 1);
    binomial[0] = 1;
    for (unsigned row = 1; row <= cells; ++row) {
        for (unsigned k = row; k != 0; --k) binomial[k] = add(binomial[k], binomial[k - 1]);
    }
    for (unsigned k = 0; k <= cells / 2; ++k) {
        full[k] = serial[k];
        full[cells - k] = (binomial[k] + MOD - serial[k]) % MOD;
    }
    for (const u32 weight : {1U, 2U}) {
        u32 value = 0;
        for (unsigned k = cells + 1; k != 0; --k) value = add(multiply(value, weight), full[k - 1]);
        require(value == evaluate_crossings(graph, weight), "full polynomial and weighted forward transfer");
        if (weight == 1) require(value == power(2, cells - 1), "half of all full boards connect each side");
    }
    std::cout << "All checks passed; supplied examples, exhaustive subsets through side 4, "
                 "game trees through side 3, pthread agreement, and weighted transfers at side 11.\n";
}

unsigned parse_unsigned(const char* text) {
    const std::string value(text);
    unsigned result = 0;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size()) {
        throw std::invalid_argument("Expected an unsigned integer");
    }
    return result;
}

}

int main(int argc, char* argv[]) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_tests();
        } else {
            if (argc > 3) throw std::invalid_argument("Usage: Euler1018 [side [threads]] | --self-test");
            const unsigned side = argc > 1 ? parse_unsigned(argv[1]) : MAX_SIDE;
            const unsigned threads = argc > 2 ? parse_unsigned(argv[2]) : 0;
            require(threads <= 64, "at most 64 threads, or 0 for automatic selection");
            const Frontier graph(side);
            std::cout << game_count(side * side, crossing_counts(graph, threads)) << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
