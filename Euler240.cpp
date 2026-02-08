#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int num_dice = 20;
    int max_points = 12;
    int top_count = 10;
    int top_sum = 70;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_int_after_prefix(arg, "--num-dice=", options.num_dice) ||
            parse_int_after_prefix(arg, "--max-points=", options.max_points) ||
            parse_int_after_prefix(arg, "--top-count=", options.top_count) ||
            parse_int_after_prefix(arg, "--top-sum=", options.top_sum)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.num_dice >= 1 && options.max_points >= 1 &&
           options.top_count >= 1 && options.top_count <= options.num_dice &&
           options.top_sum >= options.top_count && options.top_sum <= options.top_count * options.max_points;
}

std::vector<std::vector<u64>> build_binom(const int n) {
    std::vector<std::vector<u64>> c(static_cast<std::size_t>(n + 1),
                                    std::vector<u64>(static_cast<std::size_t>(n + 1), 0));
    for (int i = 0; i <= n; ++i) {
        c[static_cast<std::size_t>(i)][0] = 1;
        c[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1;
        for (int j = 1; j < i; ++j) {
            c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                c[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)] +
                c[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)];
        }
    }
    return c;
}

u64 solve(const int num_dice, const int max_points, const int top_count, const int top_sum) {
    const auto choose = build_binom(num_dice);

    std::vector<std::vector<std::vector<u64>>> dp(
        static_cast<std::size_t>(num_dice + 1),
        std::vector<std::vector<u64>>(
            static_cast<std::size_t>(top_count + 1),
            std::vector<u64>(static_cast<std::size_t>(top_sum + 1), 0)));
    dp[static_cast<std::size_t>(num_dice)][0][0] = 1;

    for (int face = max_points; face >= 1; --face) {
        std::vector<std::vector<std::vector<u64>>> next(
            static_cast<std::size_t>(num_dice + 1),
            std::vector<std::vector<u64>>(
                static_cast<std::size_t>(top_count + 1),
                std::vector<u64>(static_cast<std::size_t>(top_sum + 1), 0)));

        for (int rem = 0; rem <= num_dice; ++rem) {
            for (int filled = 0; filled <= top_count; ++filled) {
                for (int sum = 0; sum <= top_sum; ++sum) {
                    const u64 ways = dp[static_cast<std::size_t>(rem)][static_cast<std::size_t>(filled)][static_cast<std::size_t>(sum)];
                    if (ways == 0) {
                        continue;
                    }

                    for (int count = 0; count <= rem; ++count) {
                        const int add_top = std::min(count, top_count - filled);
                        const int new_filled = filled + add_top;
                        const int new_sum = sum + add_top * face;
                        if (new_sum > top_sum) {
                            continue;
                        }

                        const u64 comb = choose[static_cast<std::size_t>(rem)][static_cast<std::size_t>(count)];
                        next[static_cast<std::size_t>(rem - count)][static_cast<std::size_t>(new_filled)][static_cast<std::size_t>(new_sum)] +=
                            ways * comb;
                    }
                }
            }
        }

        dp.swap(next);
    }

    return dp[0][static_cast<std::size_t>(top_count)][static_cast<std::size_t>(top_sum)];
}

u64 brute_small(const int num_dice, const int max_points, const int top_count, const int top_sum) {
    std::vector<int> rolls(static_cast<std::size_t>(num_dice), 1);
    u64 count = 0;

    while (true) {
        std::vector<int> sorted = rolls;
        std::sort(sorted.begin(), sorted.end(), std::greater<int>());
        int sum = 0;
        for (int i = 0; i < top_count; ++i) {
            sum += sorted[static_cast<std::size_t>(i)];
        }
        if (sum == top_sum) {
            ++count;
        }

        int pos = 0;
        while (pos < num_dice && rolls[static_cast<std::size_t>(pos)] == max_points) {
            rolls[static_cast<std::size_t>(pos)] = 1;
            ++pos;
        }
        if (pos == num_dice) {
            break;
        }
        ++rolls[static_cast<std::size_t>(pos)];
    }

    return count;
}

bool run_checkpoints() {
    if (solve(5, 6, 3, 15) != 1111ULL) {
        std::cerr << "Checkpoint failed for stated sample (5,6,3,15)" << '\n';
        return false;
    }
    if (solve(7, 4, 3, 10) != brute_small(7, 4, 3, 10)) {
        std::cerr << "Checkpoint failed for brute cross-check (7,4,3,10)" << '\n';
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

    std::cout << solve(options.num_dice, options.max_points, options.top_count, options.top_sum) << '\n';
    return 0;
}
