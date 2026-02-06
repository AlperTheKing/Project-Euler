#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int kTargetIndex = 10'000'000;
constexpr u64 kSeed = 290'797ULL;
constexpr u64 kMod = 50'515'093ULL;

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

std::vector<int> generate_t_runs(int max_index) {
    std::vector<int> runs;
    runs.reserve(static_cast<std::size_t>(max_index) + 1);

    u64 s = kSeed;
    for (int i = 0; i <= max_index; ++i) {
        runs.push_back(static_cast<int>(s % 64ULL) + 1);
        s = static_cast<u64>((static_cast<u128>(s) * s) % kMod);
    }

    return runs;
}

class InvariantProfileSolver {
public:
    void process_run(int run_length, bool occupied) {
        if (run_length <= 0) {
            return;
        }

        if (occupied) {
            stack_.resize(stack_.size() + static_cast<std::size_t>(run_length), 0U);
            return;
        }

        while (run_length > 0 && !stack_.empty()) {
            close_top_pair();
            --run_length;
        }
    }

    void process_runs(const std::vector<int>& runs) {
        reset();

        bool occupied = true;
        for (const int run_length : runs) {
            process_run(run_length, occupied);
            occupied = !occupied;
        }

        finalize();
    }

    void finalize() {
        while (!stack_.empty()) {
            close_top_pair();
        }
    }

    u128 sum_of_square_final_state() const {
        u128 answer = 0;
        for (std::size_t round = 1; round < round_counts_.size(); ++round) {
            answer += static_cast<u128>(2 * round - 1) * static_cast<u128>(round_counts_[round]);
        }
        return answer;
    }

    std::vector<u64> final_state_lengths() const {
        // This is used only in small validation checkpoints.
        std::vector<u64> lengths;

        for (std::size_t len = 1; len < round_counts_.size(); ++len) {
            const u64 current = round_counts_[len];
            const u64 next = (len + 1 < round_counts_.size()) ? round_counts_[len + 1] : 0;

            if (current < next) {
                return {};
            }

            const u64 multiplicity = current - next;
            lengths.insert(lengths.end(), static_cast<std::size_t>(multiplicity),
                           static_cast<u64>(len));
        }

        return lengths;
    }

private:
    std::vector<u32> stack_;
    std::vector<u64> round_counts_{0}; // 1-indexed by elimination round.

    void reset() {
        stack_.clear();
        round_counts_.assign(1, 0);
    }

    void close_top_pair() {
        const u32 child_max_round = stack_.back();
        stack_.pop_back();

        const std::size_t round = static_cast<std::size_t>(child_max_round) + 1;
        if (round >= round_counts_.size()) {
            round_counts_.resize(round + 1, 0);
        }
        ++round_counts_[round];

        if (!stack_.empty() && stack_.back() < round) {
            stack_.back() = static_cast<u32>(round);
        }
    }
};

std::vector<int> runs_to_positions(const std::vector<int>& runs) {
    std::vector<int> positions;

    int cursor = 0;
    bool occupied = true;
    for (const int run_length : runs) {
        if (occupied) {
            positions.reserve(positions.size() + static_cast<std::size_t>(run_length));
            for (int offset = 0; offset < run_length; ++offset) {
                positions.push_back(cursor + offset);
            }
        }
        cursor += run_length;
        occupied = !occupied;
    }

    return positions;
}

std::vector<int> positions_to_runs(const std::vector<int>& positions) {
    if (positions.empty()) {
        return {};
    }

    std::vector<int> runs;

    int i = 0;
    while (i < static_cast<int>(positions.size())) {
        int j = i + 1;
        while (j < static_cast<int>(positions.size()) &&
               positions[static_cast<std::size_t>(j)] ==
                   positions[static_cast<std::size_t>(j - 1)] + 1) {
            ++j;
        }

        runs.push_back(j - i);

        if (j < static_cast<int>(positions.size())) {
            const int gap = positions[static_cast<std::size_t>(j)] -
                            positions[static_cast<std::size_t>(j - 1)] - 1;
            runs.push_back(gap);
        }

        i = j;
    }

    return runs;
}

void apply_one_turn(std::vector<int>& positions) {
    if (positions.empty()) {
        return;
    }

    std::unordered_set<int> occupied;
    occupied.reserve(positions.size() * 4 + 8);
    for (const int position : positions) {
        occupied.insert(position);
    }

    for (std::size_t i = 0; i < positions.size(); ++i) {
        const int current = positions[i];
        occupied.erase(current);

        int target = current + 1;
        while (occupied.find(target) != occupied.end()) {
            ++target;
        }

        occupied.insert(target);
        positions[i] = target;
    }
}

std::vector<int> occupied_block_lengths(const std::vector<int>& positions) {
    if (positions.empty()) {
        return {};
    }

    std::vector<int> blocks;
    int i = 0;
    while (i < static_cast<int>(positions.size())) {
        int j = i + 1;
        while (j < static_cast<int>(positions.size()) &&
               positions[static_cast<std::size_t>(j)] ==
                   positions[static_cast<std::size_t>(j - 1)] + 1) {
            ++j;
        }
        blocks.push_back(j - i);
        i = j;
    }

    return blocks;
}

std::vector<int> brute_force_final_blocks_after_turns(const std::vector<int>& runs,
                                                      int turns) {
    std::vector<int> positions = runs_to_positions(runs);
    for (int t = 0; t < turns; ++t) {
        apply_one_turn(positions);
    }
    return occupied_block_lengths(positions);
}

std::vector<u64> fast_final_blocks(const std::vector<int>& runs) {
    InvariantProfileSolver solver;
    solver.process_runs(runs);
    return solver.final_state_lengths();
}

template <typename T>
std::string vector_to_string(const std::vector<T>& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += std::to_string(values[i]);
    }
    out += "]";
    return out;
}

bool run_turn_rule_checkpoint() {
    const std::vector<int> initial = {2, 2, 2, 1, 2};
    const std::vector<int> expected = {2, 2, 1, 2, 3};

    std::vector<int> positions = runs_to_positions(initial);
    apply_one_turn(positions);
    const std::vector<int> got = positions_to_runs(positions);

    if (got != expected) {
        std::cerr << "Turn checkpoint failed: expected " << vector_to_string(expected)
                  << ", got " << vector_to_string(got) << '\n';
        return false;
    }

    return true;
}

bool run_statement_sample_checkpoint() {
    const std::vector<int> sample_runs = generate_t_runs(10);

    InvariantProfileSolver solver;
    solver.process_runs(sample_runs);

    const std::vector<u64> expected_lengths = {1, 3, 10, 24, 51, 75};
    const std::vector<u64> got_lengths = solver.final_state_lengths();

    if (got_lengths != expected_lengths) {
        std::cerr << "Sample checkpoint failed: expected final state "
                  << vector_to_string(expected_lengths) << ", got "
                  << vector_to_string(got_lengths) << '\n';
        return false;
    }

    const u64 expected_sum = 8'912;
    const u64 got_sum = static_cast<u64>(solver.sum_of_square_final_state());
    if (got_sum != expected_sum) {
        std::cerr << "Sample sum-of-squares checkpoint failed: expected " << expected_sum
                  << ", got " << got_sum << '\n';
        return false;
    }

    return true;
}

struct SmallCase {
    std::vector<int> runs;
};

std::vector<SmallCase> build_small_cross_check_cases() {
    constexpr int kCaseCount = 96;

    std::vector<SmallCase> cases;
    cases.reserve(kCaseCount);

    u64 state = 0x9e3779b97f4a7c15ULL;
    auto next_u32 = [&state]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<u32>(state >> 32);
    };

    for (int cid = 0; cid < kCaseCount; ++cid) {
        const int occupied_runs = 1 + static_cast<int>(next_u32() % 6U);

        SmallCase test_case;
        test_case.runs.reserve(static_cast<std::size_t>(2 * occupied_runs - 1));

        for (int i = 0; i < 2 * occupied_runs - 1; ++i) {
            test_case.runs.push_back(1 + static_cast<int>(next_u32() % 6U));
        }

        cases.push_back(std::move(test_case));
    }

    return cases;
}

bool run_small_cross_checks(unsigned requested_threads) {
    constexpr int kBruteTurns = 2'000;
    const std::vector<SmallCase> cases = build_small_cross_check_cases();

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = 1;
    }
    threads = std::max(1u, std::min(threads, static_cast<unsigned>(cases.size())));

    std::atomic<std::size_t> next_case{0};
    std::atomic<bool> failed{false};
    std::mutex failure_mutex;
    std::string failure_message;

    auto worker = [&]() {
        while (!failed.load(std::memory_order_relaxed)) {
            const std::size_t index = next_case.fetch_add(1, std::memory_order_relaxed);
            if (index >= cases.size()) {
                break;
            }

            const std::vector<u64> fast_u64 = fast_final_blocks(cases[index].runs);
            const std::vector<int> slow =
                brute_force_final_blocks_after_turns(cases[index].runs, kBruteTurns);

            if (fast_u64.size() != slow.size()) {
                if (!failed.exchange(true, std::memory_order_relaxed)) {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    failure_message =
                        "Random cross-check failed on case " + std::to_string(index) +
                        ": expected " + vector_to_string(slow) + ", got " +
                        vector_to_string(fast_u64) + ", runs=" +
                        vector_to_string(cases[index].runs);
                }
                break;
            }

            bool mismatch = false;
            for (std::size_t i = 0; i < slow.size(); ++i) {
                if (static_cast<u64>(slow[i]) != fast_u64[i]) {
                    mismatch = true;
                    break;
                }
            }

            if (mismatch) {
                if (!failed.exchange(true, std::memory_order_relaxed)) {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    failure_message =
                        "Random cross-check failed on case " + std::to_string(index) +
                        ": expected " + vector_to_string(slow) + ", got " +
                        vector_to_string(fast_u64) + ", runs=" +
                        vector_to_string(cases[index].runs);
                }
                break;
            }
        }
    };

    if (threads == 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve(threads);

        for (unsigned t = 0; t < threads; ++t) {
            pool.emplace_back(worker);
        }

        for (std::thread& thread : pool) {
            thread.join();
        }
    }

    if (failed.load(std::memory_order_relaxed)) {
        std::cerr << failure_message << '\n';
        return false;
    }

    return true;
}

bool run_validation_checkpoints(unsigned threads) {
    if (!run_turn_rule_checkpoint()) {
        return false;
    }
    if (!run_statement_sample_checkpoint()) {
        return false;
    }
    if (!run_small_cross_checks(threads)) {
        return false;
    }
    return true;
}

u128 solve_target() {
    InvariantProfileSolver solver;

    u64 s = kSeed;
    bool occupied = true;

    for (int index = 0; index <= kTargetIndex; ++index) {
        const int run_length = static_cast<int>(s % 64ULL) + 1;
        solver.process_run(run_length, occupied);

        occupied = !occupied;
        s = static_cast<u64>((static_cast<u128>(s) * s) % kMod);
    }

    solver.finalize();
    return solver.sum_of_square_final_state();
}

} // namespace

int main() {
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }

    if (!run_validation_checkpoints(threads)) {
        return 1;
    }

    const u128 answer = solve_target();
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
