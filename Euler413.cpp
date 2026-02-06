#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using Packed = u128;

struct PackedHash {
    std::size_t operator()(Packed value) const noexcept {
        const u64 lo = static_cast<u64>(value);
        const u64 hi = static_cast<u64>(value >> 64);
        const u64 mix = hi + 0x9e3779b97f4a7c15ULL + (lo << 6) + (lo >> 2);
        return static_cast<std::size_t>(lo ^ mix);
    }
};

using StateMap = std::unordered_map<Packed, u64, PackedHash>;

inline void decode_counts(Packed key, int d, std::array<std::uint8_t, 19>& counts) {
    for (int r = 0; r < d; ++r) {
        counts[r] = static_cast<std::uint8_t>((key >> (5 * r)) & 31U);
    }
}

inline Packed encode_counts(const std::array<std::uint8_t, 19>& counts, int d) {
    Packed key = 0;
    for (int r = 0; r < d; ++r) {
        key |= (static_cast<Packed>(counts[r]) << (5 * r));
    }
    return key;
}

inline void add_ways(StateMap& map, Packed key, u64 ways) {
    const auto it = map.find(key);
    if (it == map.end()) {
        map.emplace(key, ways);
    } else {
        it->second += ways;
    }
}

class Euler413Solver {
public:
    u64 count_for_digits(int d) const {
        // After processing a prefix, state[r] stores how many suffixes of that prefix
        // have value congruent to r (mod d). A divisible substring corresponds to r = 0.
        int next_remainder[10][19]{};
        for (int digit = 0; digit <= 9; ++digit) {
            for (int remainder = 0; remainder < d; ++remainder) {
                next_remainder[digit][remainder] = (10 * remainder + digit) % d;
            }
        }

        StateMap current_zero; // prefixes with 0 divisible substrings so far
        StateMap current_one;  // prefixes with exactly 1 divisible substring so far
        current_zero.reserve(1024);
        current_one.reserve(1024);
        current_zero.emplace(static_cast<Packed>(0), 1);

        StateMap next_zero;
        StateMap next_one;

        std::array<std::uint8_t, 19> counts{};
        std::array<std::uint8_t, 19> transformed{};

        for (int position = 1; position <= d; ++position) {
            const int min_digit = (position == 1) ? 1 : 0;

            const std::size_t expected_states =
                std::max<std::size_t>(1024,
                                      (current_zero.size() + current_one.size()) * 2 + 16);
            next_zero.clear();
            next_one.clear();
            next_zero.reserve(expected_states);
            next_one.reserve(expected_states);

            for (const auto& entry : current_zero) {
                const Packed key = entry.first;
                const u64 ways = entry.second;
                decode_counts(key, d, counts);

                for (int digit = min_digit; digit <= 9; ++digit) {
                    for (int r = 0; r < d; ++r) {
                        transformed[r] = 0;
                    }

                    for (int r = 0; r < d; ++r) {
                        const std::uint8_t count = counts[r];
                        if (count == 0) {
                            continue;
                        }
                        transformed[next_remainder[digit][r]] += count;
                    }

                    transformed[digit % d] += 1;
                    const int zero_count = transformed[0];
                    if (zero_count > 1) {
                        continue;
                    }

                    const Packed next_key = encode_counts(transformed, d);
                    if (zero_count == 0) {
                        add_ways(next_zero, next_key, ways);
                    } else {
                        add_ways(next_one, next_key, ways);
                    }
                }
            }

            for (const auto& entry : current_one) {
                const Packed key = entry.first;
                const u64 ways = entry.second;
                decode_counts(key, d, counts);

                for (int digit = min_digit; digit <= 9; ++digit) {
                    for (int r = 0; r < d; ++r) {
                        transformed[r] = 0;
                    }

                    for (int r = 0; r < d; ++r) {
                        const std::uint8_t count = counts[r];
                        if (count == 0) {
                            continue;
                        }
                        transformed[next_remainder[digit][r]] += count;
                    }

                    transformed[digit % d] += 1;
                    if (transformed[0] != 0) {
                        continue;
                    }

                    const Packed next_key = encode_counts(transformed, d);
                    add_ways(next_one, next_key, ways);
                }
            }

            current_zero.swap(next_zero);
            current_one.swap(next_one);
        }

        u64 count = 0;
        for (const auto& entry : current_one) {
            count += entry.second;
        }
        return count;
    }

    std::vector<u64> solve_all(bool allow_multithreading, unsigned requested_threads = 0) const {
        unsigned threads = 1;
        if (allow_multithreading) {
            if (requested_threads != 0) {
                threads = requested_threads;
            } else {
                const unsigned hw = std::thread::hardware_concurrency();
                threads = (hw == 0) ? 1 : hw;
            }
        }
        threads = std::max(1u, std::min(threads, 19u));

        std::vector<u64> counts(20, 0);
        if (threads == 1) {
            for (int d = 1; d <= 19; ++d) {
                counts[d] = count_for_digits(d);
            }
            return counts;
        }

        std::atomic<int> next_digit{1};
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&]() {
                while (true) {
                    const int d = next_digit.fetch_add(1, std::memory_order_relaxed);
                    if (d > 19) {
                        break;
                    }
                    counts[d] = count_for_digits(d);
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        return counts;
    }
};

bool parse_threads_argument(const std::string& arg, unsigned& threads) {
    constexpr const char* prefix = "--threads=";
    if (arg.rfind(prefix, 0) != 0) {
        return false;
    }

    const std::string value = arg.substr(std::char_traits<char>::length(prefix));
    if (value.empty()) {
        return false;
    }

    unsigned parsed = 0;
    for (const char c : value) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<unsigned>(c - '0');
    }

    threads = parsed;
    return true;
}

bool run_checkpoints(const std::vector<u64>& prefix) {
    struct Checkpoint {
        int digits;
        u64 expected;
        const char* label;
    };

    constexpr Checkpoint checks[] = {
        {1, 9, "F(10)"},
        {3, 389, "F(10^3)"},
        {7, 277674, "F(10^7)"},
    };

    for (const Checkpoint& checkpoint : checks) {
        const u64 got = prefix[checkpoint.digits];
        if (got != checkpoint.expected) {
            std::cerr << "Checkpoint failed for " << checkpoint.label
                      << ": expected " << checkpoint.expected
                      << ", got " << got << '\n';
            return false;
        }
        std::cout << checkpoint.label << " = " << got << " (ok)" << '\n';
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool allow_multithreading = true;
    unsigned requested_threads = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--single-thread") {
            allow_multithreading = false;
            continue;
        }
        if (!parse_threads_argument(arg, requested_threads)) {
            std::cerr << "Unknown argument: " << arg << '\n';
            std::cerr << "Usage: ./Euler413 [--single-thread] [--threads=N]" << '\n';
            return 1;
        }
    }

    const Euler413Solver solver;
    const std::vector<u64> by_digits = solver.solve_all(allow_multithreading, requested_threads);

    std::vector<u64> prefix(20, 0);
    for (int d = 1; d <= 19; ++d) {
        prefix[d] = prefix[d - 1] + by_digits[d];
    }

    if (!run_checkpoints(prefix)) {
        return 1;
    }

    std::cout << "F(10^19) = " << prefix[19] << '\n';
    return 0;
}
