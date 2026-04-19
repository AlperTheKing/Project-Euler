#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

constexpr u64 TARGET_N = 1'000'000'000'000'000'000ULL;
constexpr int CYCLE_START = 512;
constexpr int CYCLE_BANANAS = 71;
constexpr int CYCLE_SHIFT = 118;
constexpr int CHECK_LIMIT = 10'000;

struct Simulator {
    explicit Simulator(const int max_target)
        : offset(20 * max_target + 10'000),
          occupied(static_cast<std::size_t>(2 * offset + 10), 0) {}

    std::vector<i64> first_positions(const int max_target) {
        std::vector<i64> first(static_cast<std::size_t>(max_target + 1), 0);
        int filled = -1;

        while (filled < max_target) {
            const bool here = has(x);
            const bool next = has(x + 1);

            if (here && next) {
                clear(x + 1);
                --x;
            } else if (here) {
                clear(x);
                x += 2;
            } else if (next) {
                clear(x + 1);
                set(x);
                x += 2;
            } else {
                while (filled < std::min(line_count, max_target)) {
                    first[static_cast<std::size_t>(++filled)] = x;
                }

                assert(!has(x - 1) && !has(x) && !has(x + 1));
                set(x - 1);
                set(x);
                set(x + 1);
                x -= 2;
            }
        }

        return first;
    }

private:
    int offset;
    std::vector<unsigned char> occupied;
    i64 x = 0;
    int line_count = 0;

    std::size_t index(const i64 position) const {
        const i64 idx = position + static_cast<i64>(offset);
        assert(0 <= idx && idx < static_cast<i64>(occupied.size()));
        return static_cast<std::size_t>(idx);
    }

    bool has(const i64 position) const {
        return occupied[index(position)] != 0;
    }

    void set(const i64 position) {
        unsigned char& cell = occupied[index(position)];
        assert(cell == 0);
        cell = 1;
        ++line_count;
    }

    void clear(const i64 position) {
        unsigned char& cell = occupied[index(position)];
        assert(cell != 0);
        cell = 0;
        --line_count;
    }
};

std::string to_string_i128(i128 value) {
    if (value == 0) {
        return "0";
    }

    bool negative = false;
    if (value < 0) {
        negative = true;
        value = -value;
    }

    std::string result;
    while (value > 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) {
        result.push_back('-');
    }
    std::reverse(result.begin(), result.end());
    return result;
}

i128 bb(const u64 n, const std::vector<i64>& first) {
    const u64 target = n <= 2 ? 0 : n - 2;
    if (target < first.size()) {
        return first[static_cast<std::size_t>(target)];
    }

    const u64 cycles = (target - CYCLE_START) / CYCLE_BANANAS;
    const int residue = CYCLE_START + static_cast<int>((target - CYCLE_START) % CYCLE_BANANAS);
    return static_cast<i128>(first[static_cast<std::size_t>(residue)]) +
           static_cast<i128>(cycles) * CYCLE_SHIFT;
}

void run_checkpoints(const std::vector<i64>& first) {
    assert(bb(0, first) == 0);
    assert(bb(1, first) == 0);
    assert(bb(2, first) == 0);
    assert(bb(3, first) == 1);
    assert(bb(5, first) == -1);
    assert(bb(1000, first) == 1499);

    for (int m = CYCLE_START; m + CYCLE_BANANAS <= CHECK_LIMIT; ++m) {
        assert(first[static_cast<std::size_t>(m + CYCLE_BANANAS)] ==
               first[static_cast<std::size_t>(m)] + CYCLE_SHIFT);
    }
}

}  // namespace

int main() {
    Simulator simulator(CHECK_LIMIT);
    const std::vector<i64> first = simulator.first_positions(CHECK_LIMIT);
    run_checkpoints(first);

    std::cout << to_string_i128(bb(TARGET_N, first)) << '\n';
    return 0;
}
