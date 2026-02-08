#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    u64 target_steps = 1'000'000'000'000'000'000ULL;
    u64 detect_steps = 200'000ULL;
    u64 max_period = 500ULL;
    u64 verify_window = 50'000ULL;
    bool run_checkpoints = true;
};

struct Pattern {
    bool found = false;
    u64 start = 0ULL;
    u64 period = 0ULL;
    i64 delta = 0;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0ULL;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--target-steps=", options.target_steps) ||
            parse_u64_after_prefix(arg, "--detect-steps=", options.detect_steps) ||
            parse_u64_after_prefix(arg, "--max-period=", options.max_period) ||
            parse_u64_after_prefix(arg, "--verify-window=", options.verify_window)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.detect_steps >= 1ULL && options.max_period >= 1ULL;
}

u64 encode_cell(const std::int32_t x, const std::int32_t y) {
    return (static_cast<u64>(static_cast<std::uint32_t>(x)) << 32U) |
           static_cast<u64>(static_cast<std::uint32_t>(y));
}

std::vector<i64> simulate_black_counts(const u64 steps) {
    std::vector<i64> counts(static_cast<std::size_t>(steps + 1ULL), 0);
    std::unordered_set<u64> black;
    black.reserve(static_cast<std::size_t>(steps / 2ULL + 1024ULL));

    std::int32_t x = 0;
    std::int32_t y = 0;
    int direction = 0;  // 0=up, 1=right, 2=down, 3=left
    const int dx[4] = {0, 1, 0, -1};
    const int dy[4] = {1, 0, -1, 0};
    i64 black_count = 0;

    for (u64 step = 1ULL; step <= steps; ++step) {
        const u64 key = encode_cell(x, y);
        const auto it = black.find(key);
        if (it != black.end()) {
            black.erase(it);
            --black_count;
            direction = (direction + 3) & 3;
        } else {
            black.insert(key);
            ++black_count;
            direction = (direction + 1) & 3;
        }
        x += static_cast<std::int32_t>(dx[direction]);
        y += static_cast<std::int32_t>(dy[direction]);
        counts[static_cast<std::size_t>(step)] = black_count;
    }
    return counts;
}

Pattern detect_linear_pattern(const std::vector<i64>& counts, const u64 max_period, const u64 verify_window) {
    const u64 total_steps = static_cast<u64>(counts.size() - 1U);
    if (total_steps < 2ULL) {
        return {};
    }
    const u64 window_start = (total_steps > verify_window) ? (total_steps - verify_window) : 0ULL;

    for (u64 period = 1ULL; period <= max_period && period <= total_steps; ++period) {
        const i64 delta = counts[static_cast<std::size_t>(total_steps)] -
                          counts[static_cast<std::size_t>(total_steps - period)];
        bool tail_ok = true;
        for (u64 t = window_start; t + period <= total_steps; ++t) {
            if (counts[static_cast<std::size_t>(t + period)] -
                    counts[static_cast<std::size_t>(t)] !=
                delta) {
                tail_ok = false;
                break;
            }
        }
        if (!tail_ok) {
            continue;
        }

        const u64 upto = total_steps - period;
        std::vector<unsigned char> good(static_cast<std::size_t>(upto + 1ULL), 0U);
        for (u64 t = 0ULL; t <= upto; ++t) {
            good[static_cast<std::size_t>(t)] =
                (counts[static_cast<std::size_t>(t + period)] -
                 counts[static_cast<std::size_t>(t)] == delta)
                    ? 1U
                    : 0U;
        }

        std::vector<unsigned char> suffix_ok(static_cast<std::size_t>(upto + 2ULL), 1U);
        for (i64 t = static_cast<i64>(upto); t >= 0; --t) {
            suffix_ok[static_cast<std::size_t>(t)] =
                static_cast<unsigned char>(good[static_cast<std::size_t>(t)] &&
                                           suffix_ok[static_cast<std::size_t>(t + 1)]);
        }

        for (u64 start = 0ULL; start <= upto; ++start) {
            if (suffix_ok[static_cast<std::size_t>(start)] != 0U) {
                return Pattern{true, start, period, delta};
            }
        }
    }
    return {};
}

i64 extrapolate_black_count(const std::vector<i64>& counts, const Pattern& pattern, const u64 target_step) {
    if (target_step < counts.size()) {
        return counts[static_cast<std::size_t>(target_step)];
    }
    const u64 offset = target_step - pattern.start;
    const u64 cycles = offset / pattern.period;
    const u64 rem = offset % pattern.period;
    const i64 base = counts[static_cast<std::size_t>(pattern.start + rem)];
    const __int128 extra = static_cast<__int128>(cycles) * static_cast<__int128>(pattern.delta);
    return static_cast<i64>(static_cast<__int128>(base) + extra);
}

i64 solve(const Options& options) {
    const std::vector<i64> counts = simulate_black_counts(options.detect_steps);
    if (options.target_steps <= options.detect_steps) {
        return counts[static_cast<std::size_t>(options.target_steps)];
    }

    const Pattern pattern = detect_linear_pattern(counts, options.max_period, options.verify_window);
    if (!pattern.found) {
        throw std::runtime_error("Could not detect linear periodic tail");
    }
    return extrapolate_black_count(counts, pattern, options.target_steps);
}

bool run_checkpoints() {
    const std::vector<i64> small = simulate_black_counts(200ULL);
    if (small[10] != 6) {
        std::cerr << "Checkpoint failed for first 10 moves" << '\n';
        return false;
    }

    const std::vector<i64> large = simulate_black_counts(200'000ULL);
    const Pattern pattern = detect_linear_pattern(large, 500ULL, 50'000ULL);
    if (!pattern.found || pattern.period != 104ULL || pattern.delta != 12LL) {
        std::cerr << "Checkpoint failed for detected highway period/delta" << '\n';
        return false;
    }

    const i64 predicted = extrapolate_black_count(large, pattern, 150'000ULL);
    if (predicted != large[150'000]) {
        std::cerr << "Checkpoint failed for extrapolation consistency at step 150000" << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    try {
        std::cout << solve(options) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }
    return 0;
}
