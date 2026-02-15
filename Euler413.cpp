#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;

using StateMap = std::unordered_map<u64, u64>;

u64 encode_state(const std::array<std::uint8_t, 19>& counts, int n) {
    u64 key = 0;
    u64 mul = 1;
    for (int i = 0; i < n; ++i) {
        key += static_cast<u64>(counts[i]) * mul;
        mul *= 3ULL;
    }
    return key;
}

void decode_state(u64 key, int n, std::array<std::uint8_t, 19>& counts) {
    for (int i = 0; i < n; ++i) {
        counts[i] = static_cast<std::uint8_t>(key % 3ULL);
        key /= 3ULL;
    }
}

inline void add_ways(StateMap& map, u64 key, u64 ways) {
    const auto it = map.find(key);
    if (it == map.end()) {
        map.emplace(key, ways);
    } else {
        it->second += ways;
    }
}

u64 count_for_length(int n) {
    std::array<std::array<int, 19>, 10> next_rem{};
    std::array<int, 10> digit_mod{};
    for (int d = 0; d <= 9; ++d) {
        digit_mod[d] = d % n;
        for (int r = 0; r < n; ++r) {
            next_rem[d][r] = (10 * r + digit_mod[d]) % n;
        }
    }

    StateMap curr0;
    StateMap curr1;
    curr0.reserve(1024);
    curr1.reserve(1024);

    std::array<std::uint8_t, 19> counts{};
    for (int d = 1; d <= 9; ++d) {
        std::fill(counts.begin(), counts.begin() + n, 0);
        counts[digit_mod[d]] = 1;
        const u64 key = encode_state(counts, n);
        if (counts[0] == 0) {
            add_ways(curr0, key, 1);
        } else if (counts[0] == 1) {
            add_ways(curr1, key, 1);
        }
    }

    StateMap next0;
    StateMap next1;

    std::array<std::uint8_t, 19> next_counts{};

    for (int pos = 2; pos <= n; ++pos) {
        const std::size_t expected_states =
            std::max<std::size_t>(1024, (curr0.size() + curr1.size()) * 2 + 16);
        next0.clear();
        next1.clear();
        next0.reserve(expected_states);
        next1.reserve(expected_states);

        for (const auto& entry : curr0) {
            const u64 key = entry.first;
            const u64 ways = entry.second;
            decode_state(key, n, counts);

            for (int d = 0; d <= 9; ++d) {
                std::fill(next_counts.begin(), next_counts.begin() + n, 0);
                for (int r = 0; r < n; ++r) {
                    const std::uint8_t c = counts[r];
                    if (!c) continue;
                    const int nr = next_rem[d][r];
                    const std::uint8_t sum = static_cast<std::uint8_t>(
                        std::min<int>(2, next_counts[nr] + c));
                    next_counts[nr] = sum;
                }
                const int dr = digit_mod[d];
                next_counts[dr] = static_cast<std::uint8_t>(std::min<int>(2, next_counts[dr] + 1));

                const int zeros = next_counts[0];
                if (zeros == 0) {
                    add_ways(next0, encode_state(next_counts, n), ways);
                } else if (zeros == 1) {
                    add_ways(next1, encode_state(next_counts, n), ways);
                }
            }
        }

        for (const auto& entry : curr1) {
            const u64 key = entry.first;
            const u64 ways = entry.second;
            decode_state(key, n, counts);

            for (int d = 0; d <= 9; ++d) {
                std::fill(next_counts.begin(), next_counts.begin() + n, 0);
                for (int r = 0; r < n; ++r) {
                    const std::uint8_t c = counts[r];
                    if (!c) continue;
                    const int nr = next_rem[d][r];
                    const std::uint8_t sum = static_cast<std::uint8_t>(
                        std::min<int>(2, next_counts[nr] + c));
                    next_counts[nr] = sum;
                }
                const int dr = digit_mod[d];
                next_counts[dr] = static_cast<std::uint8_t>(std::min<int>(2, next_counts[dr] + 1));

                if (next_counts[0] == 0) {
                    add_ways(next1, encode_state(next_counts, n), ways);
                }
            }
        }

        curr0.swap(next0);
        curr1.swap(next1);
    }

    u64 total = 0;
    for (const auto& entry : curr1) {
        total += entry.second;
    }
    return total;
}

bool run_checkpoints(const std::vector<u64>& prefix) {
    struct Checkpoint {
        int digits;
        u64 expected;
    };
    const Checkpoint checks[] = {
        {1, 9},
        {3, 389},
        {7, 277674},
    };
    for (const auto& checkpoint : checks) {
        if (prefix[checkpoint.digits] != checkpoint.expected) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    std::vector<u64> counts(20, 0);
    for (int n = 1; n <= 19; ++n) {
        counts[n] = count_for_length(n);
    }

    std::vector<u64> prefix(20, 0);
    for (int n = 1; n <= 19; ++n) {
        prefix[n] = prefix[n - 1] + counts[n];
    }

    if (!run_checkpoints(prefix)) {
        std::cerr << "Checkpoint failed\n";
        return 1;
    }

    std::cout << "F(10^19) = " << prefix[19] << '\n';
    return 0;
}
