#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u8 = std::uint8_t;
using u64 = std::uint64_t;

struct Options {
    int limit = 1000;
    bool run_checkpoints = true;
    int checkpoint_bruteforce_limit = 35;
    unsigned requested_threads = 0U;
};

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (static_cast<u64>(std::numeric_limits<unsigned>::max()) - digit) /
                          10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_int_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    std::int64_t parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<std::int64_t>(c - '0');
        if (parsed > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
            return false;
        }
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--check-limit=",
                                   options.checkpoint_bruteforce_limit)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

std::size_t checked_square_size(const int side) {
    if (side <= 0) {
        throw std::runtime_error("Side length must be positive");
    }

    const std::size_t s = static_cast<std::size_t>(side);
    if (s > std::numeric_limits<std::size_t>::max() / s) {
        throw std::overflow_error("2D table size overflow");
    }
    return s * s;
}

std::size_t checked_cube_size(const int side) {
    const std::size_t square = checked_square_size(side);
    const std::size_t s = static_cast<std::size_t>(side);
    if (square > std::numeric_limits<std::size_t>::max() / s) {
        throw std::overflow_error("3D table size overflow");
    }
    return square * s;
}

u64 solve_fast(const int limit) {
    if (limit < 0) {
        throw std::invalid_argument("Limit must be non-negative");
    }

    const int side = limit + 1;
    const std::size_t table_size = checked_square_size(side);

    std::vector<u8> one(table_size, 0U);
    std::vector<u8> two(table_size, 0U);
    std::vector<u8> all(table_size, 0U);

    const auto id = [side](const int a, const int b) -> std::size_t {
        return static_cast<std::size_t>(a) * static_cast<std::size_t>(side) +
               static_cast<std::size_t>(b);
    };

    u64 sum = 0ULL;

    for (int x = 0; x <= limit; ++x) {
        for (int y = x; y <= limit; ++y) {
            for (int z = y; z <= limit; ++z) {
                if (one[id(y, z)] != 0U || one[id(x, z)] != 0U ||
                    one[id(x, y)] != 0U || two[id(y - x, z)] != 0U ||
                    two[id(z - y, x)] != 0U || two[id(z - x, y)] != 0U ||
                    all[id(y - x, z - x)] != 0U) {
                    continue;
                }

                sum += static_cast<u64>(x + y + z);

                one[id(y, z)] = 1U;
                one[id(x, z)] = 1U;
                one[id(x, y)] = 1U;

                two[id(y - x, z)] = 1U;
                two[id(z - y, x)] = 1U;
                two[id(z - x, y)] = 1U;

                all[id(y - x, z - x)] = 1U;
            }
        }
    }

    return sum;
}

inline void sort_triplet(int& a, int& b, int& c) {
    if (a > b) {
        std::swap(a, b);
    }
    if (b > c) {
        std::swap(b, c);
    }
    if (a > b) {
        std::swap(a, b);
    }
}

struct BruteResult {
    int limit = 0;
    int side = 0;
    u64 losing_sum = 0ULL;
    std::vector<u8> is_losing;

    bool query_is_losing(int x, int y, int z) const {
        sort_triplet(x, y, z);
        const std::size_t index =
            (static_cast<std::size_t>(x) * static_cast<std::size_t>(side) +
             static_cast<std::size_t>(y)) *
                static_cast<std::size_t>(side) +
            static_cast<std::size_t>(z);
        return is_losing[index] != 0U;
    }
};

bool has_move_to_losing_state(const std::vector<u8>& is_losing,
                              const int side,
                              const int x,
                              const int y,
                              const int z) {
    const auto index3 = [side](const int a, const int b, const int c) -> std::size_t {
        return (static_cast<std::size_t>(a) * static_cast<std::size_t>(side) +
                static_cast<std::size_t>(b)) *
                   static_cast<std::size_t>(side) +
               static_cast<std::size_t>(c);
    };

    const auto can_reach_losing = [&is_losing, &index3](int a, int b, int c) {
        sort_triplet(a, b, c);
        return is_losing[index3(a, b, c)] != 0U;
    };

    for (int n = 1; n <= x; ++n) {
        if (can_reach_losing(x - n, y, z)) {
            return true;
        }
    }
    for (int n = 1; n <= y; ++n) {
        if (can_reach_losing(x, y - n, z)) {
            return true;
        }
    }
    for (int n = 1; n <= z; ++n) {
        if (can_reach_losing(x, y, z - n)) {
            return true;
        }
    }

    for (int n = 1; n <= x; ++n) {
        if (can_reach_losing(x - n, y - n, z)) {
            return true;
        }
    }
    for (int n = 1; n <= x; ++n) {
        if (can_reach_losing(x - n, y, z - n)) {
            return true;
        }
    }
    for (int n = 1; n <= y; ++n) {
        if (can_reach_losing(x, y - n, z - n)) {
            return true;
        }
    }

    for (int n = 1; n <= x; ++n) {
        if (can_reach_losing(x - n, y - n, z - n)) {
            return true;
        }
    }

    return false;
}

BruteResult solve_bruteforce_parallel(const int limit, unsigned thread_count) {
    if (limit < 0) {
        throw std::invalid_argument("Brute-force limit must be non-negative");
    }

    const int side = limit + 1;
    BruteResult result;
    result.limit = limit;
    result.side = side;
    result.is_losing.assign(checked_cube_size(side), 0U);

    const auto index3 = [side](const int a, const int b, const int c) -> std::size_t {
        return (static_cast<std::size_t>(a) * static_cast<std::size_t>(side) +
                static_cast<std::size_t>(b)) *
                   static_cast<std::size_t>(side) +
               static_cast<std::size_t>(c);
    };

    thread_count = std::max(1U, thread_count);

    for (int total = 0; total <= 3 * limit; ++total) {
        std::vector<std::array<int, 3>> layer_states;

        const int x_hi = std::min(limit, total);
        for (int x = 0; x <= x_hi; ++x) {
            const int rem = total - x;
            const int y_hi = std::min(limit, rem / 2);
            for (int y = x; y <= y_hi; ++y) {
                const int z = rem - y;
                if (z < y || z > limit) {
                    continue;
                }
                layer_states.push_back({x, y, z});
            }
        }

        if (layer_states.empty()) {
            continue;
        }

        const unsigned workers =
            std::min<unsigned>(thread_count, static_cast<unsigned>(layer_states.size()));

        struct ChunkResult {
            u64 partial_sum = 0ULL;
            std::vector<std::size_t> losing_indices;
        };

        std::vector<ChunkResult> chunks(workers);

        const std::size_t per_worker =
            (layer_states.size() + static_cast<std::size_t>(workers) - 1U) /
            static_cast<std::size_t>(workers);

        auto worker_fn = [&](const unsigned worker_id) {
            const std::size_t begin = static_cast<std::size_t>(worker_id) * per_worker;
            const std::size_t end =
                std::min(layer_states.size(), begin + per_worker);

            ChunkResult& chunk = chunks[worker_id];
            chunk.losing_indices.reserve(end > begin ? (end - begin) / 4 + 1 : 1);

            for (std::size_t i = begin; i < end; ++i) {
                const int x = layer_states[i][0];
                const int y = layer_states[i][1];
                const int z = layer_states[i][2];

                if (has_move_to_losing_state(result.is_losing, side, x, y, z)) {
                    continue;
                }

                chunk.partial_sum += static_cast<u64>(x + y + z);
                chunk.losing_indices.push_back(index3(x, y, z));
            }
        };

        std::vector<std::thread> pool;
        pool.reserve(workers > 0U ? workers - 1U : 0U);

        for (unsigned worker = 1U; worker < workers; ++worker) {
            pool.emplace_back(worker_fn, worker);
        }

        worker_fn(0U);

        for (auto& t : pool) {
            t.join();
        }

        for (const ChunkResult& chunk : chunks) {
            result.losing_sum += chunk.partial_sum;
            for (const std::size_t idx : chunk.losing_indices) {
                result.is_losing[idx] = 1U;
            }
        }
    }

    return result;
}

bool check_equal_u64(const u64 actual, const u64 expected, const std::string& label) {
    if (actual == expected) {
        return true;
    }

    std::cerr << "Checkpoint failed for " << label << ": expected " << expected
              << ", got " << actual << '\n';
    return false;
}

bool check_bool(const bool actual, const bool expected, const std::string& label) {
    if (actual == expected) {
        return true;
    }

    std::cerr << "Checkpoint failed for " << label << ": expected "
              << (expected ? "true" : "false") << ", got "
              << (actual ? "true" : "false") << '\n';
    return false;
}

bool run_checkpoints(const unsigned thread_count, const int brute_limit) {
    if (brute_limit < 13) {
        std::cerr << "Checkpoint brute-force limit must be at least 13.\n";
        return false;
    }

    const BruteResult brute = solve_bruteforce_parallel(brute_limit, thread_count);
    const u64 fast_at_brute_limit = solve_fast(brute_limit);
    if (!check_equal_u64(fast_at_brute_limit, brute.losing_sum,
                         "fast vs brute at n=" + std::to_string(brute_limit))) {
        return false;
    }

    if (!check_bool(brute.query_is_losing(0, 1, 2), true,
                    "losing example (0,1,2)")) {
        return false;
    }
    if (!check_bool(brute.query_is_losing(1, 3, 3), true,
                    "losing example (1,3,3)")) {
        return false;
    }
    if (!check_bool(brute.query_is_losing(0, 0, 13), false,
                    "winning example (0,0,13)")) {
        return false;
    }
    if (!check_bool(brute.query_is_losing(0, 11, 11), false,
                    "winning example (0,11,11)")) {
        return false;
    }
    if (!check_bool(brute.query_is_losing(5, 5, 5), false,
                    "winning example (5,5,5)")) {
        return false;
    }

    const u64 sample = solve_fast(100);
    if (!check_equal_u64(sample, 173895ULL, "sample sum at n=100")) {
        return false;
    }

    return true;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }

        if (options.limit < 0) {
            std::cerr << "--limit must be non-negative.\n";
            return 1;
        }
        if (options.checkpoint_bruteforce_limit < 0) {
            std::cerr << "--check-limit must be non-negative.\n";
            return 1;
        }

        const unsigned thread_count = pick_thread_count(options.requested_threads);

        if (options.run_checkpoints &&
            !run_checkpoints(thread_count, options.checkpoint_bruteforce_limit)) {
            return 1;
        }

        const u64 answer = solve_fast(options.limit);
        std::cout << answer << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
