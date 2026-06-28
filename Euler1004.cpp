#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = std::int64_t;

constexpr i64 MOD = 1'000'000'007LL;
constexpr int DIGITS = 10;
constexpr int MAX_CELLS = DIGITS * DIGITS;

std::array<i64, MAX_CELLS + 1> factorial;

i64 mod_pow(i64 base, i64 exp) {
    i64 result = 1;
    while (exp > 0) {
        if ((exp & 1LL) != 0) {
            result = result * base % MOD;
        }
        base = base * base % MOD;
        exp >>= 1LL;
    }
    return result;
}

i64 mod_inverse(const i64 value) {
    return mod_pow(value, MOD - 2);
}

i64 shape_word_count(const std::vector<int>& partition) {
    int cells = 0;
    for (const int row : partition) {
        cells += row;
    }

    i64 hook_product = 1;
    i64 content_product = 1;
    for (int i = 0; i < static_cast<int>(partition.size()); ++i) {
        for (int j = 0; j < partition[i]; ++j) {
            int below = 0;
            for (int r = i + 1; r < static_cast<int>(partition.size()); ++r) {
                if (partition[r] > j) {
                    ++below;
                }
            }

            const int right = partition[i] - j - 1;
            const int hook = right + below + 1;
            const int content = DIGITS + (j + 1) - (i + 1);
            hook_product = hook_product * hook % MOD;
            content_product = content_product * content % MOD;
        }
    }

    const i64 inv_hooks = mod_inverse(hook_product);
    return factorial[cells] * content_product % MOD * inv_hooks % MOD * inv_hooks % MOD;
}

void enumerate_partitions(const int max_part,
                          const int max_rows,
                          std::vector<int>& partition,
                          const int max_cells,
                          i64& balanced,
                          i64& decreasing_excess) {
    if (!partition.empty()) {
        int cells = 0;
        for (const int row : partition) {
            cells += row;
        }

        if (cells <= max_cells) {
            const i64 ways = shape_word_count(partition);
            const int width = partition.front();
            const int height = static_cast<int>(partition.size());
            if (width == height) {
                balanced += ways;
                if (balanced >= MOD) {
                    balanced -= MOD;
                }
            }
            if (height == width + 1) {
                decreasing_excess += ways;
                if (decreasing_excess >= MOD) {
                    decreasing_excess -= MOD;
                }
            }
        }
    }

    if (static_cast<int>(partition.size()) == max_rows) {
        return;
    }

    int used = 0;
    for (const int row : partition) {
        used += row;
    }

    for (int next = max_part; next >= 1; --next) {
        if (used + next > max_cells) {
            continue;
        }
        partition.push_back(next);
        enumerate_partitions(next, max_rows, partition, max_cells, balanced, decreasing_excess);
        partition.pop_back();
    }
}

std::pair<i64, i64> count_all_words(const int max_cells) {
    i64 balanced = 0;
    i64 decreasing_excess = 0;
    std::vector<int> partition;
    enumerate_partitions(DIGITS, DIGITS, partition, max_cells, balanced, decreasing_excess);
    return {balanced, decreasing_excess};
}

i64 positive_balanced_count(const int max_digits) {
    const auto [balanced, ignored] = count_all_words(max_digits);
    const auto [unused, decreasing_excess] = count_all_words(max_digits - 1);
    (void)ignored;
    (void)unused;

    i64 result = (balanced - decreasing_excess - 1) % MOD;
    if (result < 0) {
        result += MOD;
    }
    return result;
}

void run_checkpoints() {
    assert(positive_balanced_count(4) == 2274);
}

}  // namespace

int main() {
    factorial[0] = 1;
    for (int i = 1; i <= MAX_CELLS; ++i) {
        factorial[i] = factorial[i - 1] * i % MOD;
    }

    run_checkpoints();
    std::cout << positive_balanced_count(MAX_CELLS) << '\n';
    return 0;
}
