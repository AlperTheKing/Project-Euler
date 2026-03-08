#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int STARTS = 10;
constexpr int RANKS = 13;
constexpr int MAX_STRAIGHTS = 8;
constexpr int FREE_STATE = 4;

using Pattern = std::array<int, STARTS>;

constexpr std::array<std::array<int, 5>, STARTS> COVER{{
    {{0, 1, 2, 3, 12}},
    {{0, 1, 2, 3, 4}},
    {{1, 2, 3, 4, 5}},
    {{2, 3, 4, 5, 6}},
    {{3, 4, 5, 6, 7}},
    {{4, 5, 6, 7, 8}},
    {{5, 6, 7, 8, 9}},
    {{6, 7, 8, 9, 10}},
    {{7, 8, 9, 10, 11}},
    {{8, 9, 10, 11, 12}},
}};

constexpr std::array<u64, MAX_STRAIGHTS + 1> FACT{{
    1ULL,
    1ULL,
    2ULL,
    6ULL,
    24ULL,
    120ULL,
    720ULL,
    5040ULL,
    40320ULL,
}};

constexpr std::array<int, MAX_STRAIGHTS + 1> POW5{{
    1,
    5,
    25,
    125,
    625,
    3125,
    15625,
    78125,
    390625,
}};

const std::array<std::vector<std::array<int, 4>>, 5> SUIT_PERMS = [] {
    std::array<std::vector<std::array<int, 4>>, 5> perms;
    perms[0].push_back({0, 0, 0, 0});
    for (int a = 0; a < 4; ++a) {
        perms[1].push_back({a, 0, 0, 0});
        for (int b = 0; b < 4; ++b) {
            if (b == a) continue;
            perms[2].push_back({a, b, 0, 0});
            for (int c = 0; c < 4; ++c) {
                if (c == a || c == b) continue;
                perms[3].push_back({a, b, c, 0});
                for (int d = 0; d < 4; ++d) {
                    if (d == a || d == b || d == c) continue;
                    perms[4].push_back({a, b, c, d});
                }
            }
        }
    }
    return perms;
}();

struct PatternData {
    int straight_count = 0;
    std::array<int, MAX_STRAIGHTS> start{};
    std::array<int, MAX_STRAIGHTS> first_rank{};
    std::array<int, MAX_STRAIGHTS> last_rank{};
    std::array<std::array<bool, RANKS>, MAX_STRAIGHTS> active{};
    std::array<std::array<int, MAX_STRAIGHTS>, RANKS> entry_pos{};
    std::array<std::array<int, MAX_STRAIGHTS>, RANKS> next_pos{};
    std::array<std::array<int, MAX_STRAIGHTS>, RANKS> current_labels{};
    std::array<std::array<int, MAX_STRAIGHTS>, RANKS> next_labels{};
    std::array<int, RANKS> current_count{};
    std::array<int, RANKS> next_count{};
    u64 divisor = 1;
};

PatternData build_pattern_data(const Pattern& pattern) {
    PatternData data;
    for (int s = 0; s < STARTS; ++s) {
        for (int k = 0; k < pattern[s]; ++k) {
            const int label = data.straight_count++;
            data.start[label] = s;
            data.first_rank[label] = COVER[static_cast<std::size_t>(s)][0];
            data.last_rank[label] = COVER[static_cast<std::size_t>(s)][4];
            for (const int rank : COVER[static_cast<std::size_t>(s)]) {
                data.active[label][rank] = true;
            }
        }
        data.divisor *= FACT[static_cast<std::size_t>(pattern[s])];
    }

    for (int rank = 0; rank < RANKS; ++rank) {
        data.entry_pos[rank].fill(-1);
        data.next_pos[rank].fill(-1);

        int current_idx = 0;
        int next_idx = 0;
        for (int label = 0; label < data.straight_count; ++label) {
            if (data.active[label][rank]) {
                data.current_labels[rank][current_idx++] = label;
            }
            if (data.first_rank[label] <= rank && rank < data.last_rank[label]) {
                data.next_pos[rank][label] = next_idx;
                data.next_labels[rank][next_idx++] = label;
            }
        }
        data.current_count[rank] = current_idx;
        data.next_count[rank] = next_idx;

        int entry_idx = 0;
        for (int label = 0; label < data.straight_count; ++label) {
            if (data.first_rank[label] < rank && rank <= data.last_rank[label]) {
                data.entry_pos[rank][label] = entry_idx++;
            }
        }
    }

    return data;
}

u64 count_pattern(const Pattern& pattern) {
    const PatternData data = build_pattern_data(pattern);

    std::unordered_map<int, u64> current;
    std::unordered_map<int, u64> next;
    current.reserve(1024);
    next.reserve(1024);
    current.emplace(0, 1ULL);

    for (int rank = 0; rank < RANKS; ++rank) {
        next.clear();

        const int current_count = data.current_count[rank];
        const int next_count = data.next_count[rank];

        for (const auto& [state_code, ways] : current) {
            std::array<int, MAX_STRAIGHTS> state_digits{};
            int tmp = state_code;
            for (int i = 0; i < MAX_STRAIGHTS; ++i) {
                state_digits[i] = tmp % 5;
                tmp /= 5;
            }

            u64 carried_code = 0;
            for (int idx = 0; idx < next_count; ++idx) {
                const int label = data.next_labels[rank][idx];
                if (!data.active[label][rank]) {
                    const int old_pos = data.entry_pos[rank][label];
                    carried_code += static_cast<u64>(state_digits[old_pos]) *
                                    static_cast<u64>(POW5[idx]);
                }
            }

            for (const auto& suits : SUIT_PERMS[static_cast<std::size_t>(current_count)]) {
                bool ok = true;
                u64 next_code = carried_code;

                for (int idx = 0; idx < current_count; ++idx) {
                    const int label = data.current_labels[rank][idx];
                    const int suit = suits[static_cast<std::size_t>(idx)];
                    const int old_pos = data.entry_pos[rank][label];

                    int new_state = suit;
                    if (old_pos != -1) {
                        const int old_state = state_digits[old_pos];
                        new_state = (old_state == FREE_STATE || old_state != suit) ? FREE_STATE : old_state;
                    }

                    if (rank == data.last_rank[label]) {
                        if (new_state != FREE_STATE) {
                            ok = false;
                            break;
                        }
                    } else {
                        const int out_pos = data.next_pos[rank][label];
                        next_code += static_cast<u64>(new_state) *
                                     static_cast<u64>(POW5[out_pos]);
                    }
                }

                if (!ok) continue;
                auto it = next.find(static_cast<int>(next_code));
                if (it == next.end()) {
                    next.emplace(static_cast<int>(next_code), ways);
                } else {
                    it->second += ways;
                }
            }
        }

        current.swap(next);
    }

    const auto it = current.find(0);
    assert(it != current.end());
    return it->second / data.divisor;
}

void generate_patterns_rec(int start_idx,
                           int remaining,
                           Pattern& pattern,
                           std::array<int, RANKS>& occupancy,
                           std::vector<Pattern>& out) {
    if (start_idx == STARTS) {
        if (remaining == 0) out.push_back(pattern);
        return;
    }

    const auto& ranks = COVER[static_cast<std::size_t>(start_idx)];
    for (int copies = 0; copies <= std::min(4, remaining); ++copies) {
        bool ok = true;
        for (const int rank : ranks) {
            if (occupancy[rank] + copies > 4) {
                ok = false;
                break;
            }
        }
        if (!ok) break;

        pattern[start_idx] = copies;
        for (const int rank : ranks) occupancy[rank] += copies;
        generate_patterns_rec(start_idx + 1, remaining - copies, pattern, occupancy, out);
        for (const int rank : ranks) occupancy[rank] -= copies;
    }
    pattern[start_idx] = 0;
}

std::vector<Pattern> generate_patterns(int straight_count) {
    std::vector<Pattern> patterns;
    Pattern pattern{};
    std::array<int, RANKS> occupancy{};
    generate_patterns_rec(0, straight_count, pattern, occupancy, patterns);
    return patterns;
}

u64 count_total_single_thread(int straight_count) {
    const std::vector<Pattern> patterns = generate_patterns(straight_count);
    u128 total = 0;
    for (const Pattern& pattern : patterns) {
        total += count_pattern(pattern);
    }
    return static_cast<u64>(total);
}

struct WorkerCtx {
    const std::vector<Pattern>* patterns = nullptr;
    std::atomic<std::size_t>* next_index = nullptr;
    u128 subtotal = 0;
};

void* worker_main(void* raw_ctx) {
    auto* ctx = static_cast<WorkerCtx*>(raw_ctx);
    for (;;) {
        const std::size_t idx = ctx->next_index->fetch_add(1, std::memory_order_relaxed);
        if (idx >= ctx->patterns->size()) break;
        ctx->subtotal += count_pattern((*ctx->patterns)[idx]);
    }
    return nullptr;
}

u64 count_total_parallel(int straight_count) {
    const std::vector<Pattern> patterns = generate_patterns(straight_count);

    long cpu_count = ::sysconf(_SC_NPROCESSORS_ONLN);
    if (cpu_count < 1) cpu_count = 1;
    const int thread_count = static_cast<int>(std::min<long>(32, cpu_count));

    std::atomic<std::size_t> next_index(0);
    std::vector<pthread_t> threads(static_cast<std::size_t>(thread_count));
    std::vector<WorkerCtx> workers(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        workers[static_cast<std::size_t>(t)].patterns = &patterns;
        workers[static_cast<std::size_t>(t)].next_index = &next_index;
        const int rc = pthread_create(&threads[static_cast<std::size_t>(t)], nullptr, worker_main,
                                      &workers[static_cast<std::size_t>(t)]);
        assert(rc == 0);
    }
    for (int t = 0; t < thread_count; ++t) {
        const int rc = pthread_join(threads[static_cast<std::size_t>(t)], nullptr);
        assert(rc == 0);
    }

    u128 total = 0;
    for (const WorkerCtx& worker : workers) {
        total += worker.subtotal;
    }
    return static_cast<u64>(total);
}

std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    assert(count_total_single_thread(1) == 10200ULL);
    assert(count_total_single_thread(2) == 31832952ULL);

    const u64 answer = count_total_parallel(8);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
