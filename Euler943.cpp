#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kN = 22'332'223'332'233ULL;
constexpr u64 kMod = 2'233'222'333ULL;
constexpr int kMinValue = 2;
constexpr int kMaxValue = 223;
constexpr int kMaxDepth = 61;

struct Summary {
    u64 count_a = 0;
    u64 count_b = 0;
    u64 next_state = 0;

    u64 total() const {
        return count_a + count_b;
    }
};

class KolakoskiPrefixCounter {
public:
    KolakoskiPrefixCounter(int a, int b) : a_(a), b_(b) {
        memo_.reserve(1U << 15);
    }

    Summary evaluate_prefix(u64 n) {
        for (int depth = 0; depth <= kMaxDepth; ++depth) {
            Summary s = evaluate_block(0, depth, n);
            if (s.total() >= n) {
                return s;
            }
        }
        assert(false);
        return {};
    }

private:
    int a_;
    int b_;
    std::unordered_map<u64, Summary> memo_;

    static u64 cache_key(u64 state, int depth) {
        return state + (4ULL << depth);
    }

    Summary evaluate_block(u64 state, int depth, u64 budget) {
        const u64 length_mask = 2ULL << depth;
        const u64 bit = state & length_mask;
        const u64 run_length = bit ? static_cast<u64>(b_) : static_cast<u64>(a_);
        const u64 take = std::min(run_length, budget);

        if (depth == 0) {
            Summary out;
            if ((state & 1ULL) == 0ULL) {
                out.count_a = take;
            } else {
                out.count_b = take;
            }
            out.next_state = state ^ 1ULL;
            return out;
        }

        u64 produced_a = 0;
        u64 produced_b = 0;
        u64 child_state = state ^ bit;

        for (u64 i = 0; i < take; ++i) {
            const u64 key = cache_key(child_state, depth - 1);
            Summary child;
            const auto it = memo_.find(key);
            if (it != memo_.end()) {
                const u64 cached_total = it->second.total();
                if (produced_a + produced_b + cached_total <= budget) {
                    child = it->second;
                } else {
                    child = evaluate_block(child_state, depth - 1, budget - (produced_a + produced_b));
                }
            } else {
                child = evaluate_block(child_state, depth - 1, budget - (produced_a + produced_b));
            }
            produced_a += child.count_a;
            produced_b += child.count_b;
            child_state = child.next_state;
        }

        Summary out;
        out.count_a = produced_a;
        out.count_b = produced_b;
        out.next_state = child_state ^ bit ^ (1ULL << depth);
        memo_[cache_key(state, depth)] = out;
        return out;
    }
};

u64 compute_t(int a, int b, u64 n) {
    KolakoskiPrefixCounter counter(a, b);
    const Summary s = counter.evaluate_prefix(n);
    const u128 total = static_cast<u128>(s.count_a) * static_cast<u64>(a) +
                       static_cast<u128>(s.count_b) * static_cast<u64>(b);
    return static_cast<u64>(total);
}

u64 brute_t(int a, int b, u64 n) {
    std::vector<int> seq;
    seq.reserve(static_cast<std::size_t>(n + static_cast<u64>(b) + 8));
    seq.push_back(a);

    std::size_t run_idx = 0;
    int current_symbol = a;
    int remaining = a - 1;
    u64 count_a = 1;
    u64 count_b = 0;

    while (seq.size() < n) {
        if (remaining == 0) {
            ++run_idx;
            current_symbol = (run_idx % 2 == 0) ? a : b;
            remaining = seq[run_idx];
        }
        seq.push_back(current_symbol);
        if (current_symbol == a) {
            ++count_a;
        } else {
            ++count_b;
        }
        --remaining;
    }

    return count_a * static_cast<u64>(a) + count_b * static_cast<u64>(b);
}

void run_validations() {
    assert(compute_t(2, 3, 10) == 25ULL);
    assert(compute_t(3, 4, 20) == 70ULL);

    for (int a = 2; a <= 7; ++a) {
        for (int b = 2; b <= 7; ++b) {
            if (a == b) {
                continue;
            }
            for (u64 n = 1; n <= 180; n += 11) {
                assert(compute_t(a, b, n) == brute_t(a, b, n));
            }
        }
    }
}

u64 solve_parallel() {
    const int count_a_values = kMaxValue - kMinValue + 1;
    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) {
        threads = 4;
    }
    threads = std::min<unsigned int>(threads, static_cast<unsigned int>(count_a_values));

    std::atomic<int> next_a{kMinValue};
    std::vector<u64> partial(threads, 0);
    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned int tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            u64 local = 0;
            while (true) {
                const int a = next_a.fetch_add(1);
                if (a > kMaxValue) {
                    break;
                }
                for (int b = kMinValue; b <= kMaxValue; ++b) {
                    if (a == b) {
                        continue;
                    }
                    local += compute_t(a, b, kN) % kMod;
                    if (local >= kMod) {
                        local %= kMod;
                    }
                }
            }
            partial[tid] = local % kMod;
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    u64 answer = 0;
    for (u64 part : partial) {
        answer += part;
        answer %= kMod;
    }
    return answer;
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve_parallel() << '\n';
    return 0;
}
