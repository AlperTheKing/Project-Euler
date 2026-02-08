#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int pete_dice = 9;
    int pete_sides = 4;
    int colin_dice = 6;
    int colin_sides = 6;
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

        if (parse_int_after_prefix(arg, "--pete-dice=", options.pete_dice) ||
            parse_int_after_prefix(arg, "--pete-sides=", options.pete_sides) ||
            parse_int_after_prefix(arg, "--colin-dice=", options.colin_dice) ||
            parse_int_after_prefix(arg, "--colin-sides=", options.colin_sides)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.pete_dice >= 1 && options.colin_dice >= 1 &&
           options.pete_sides >= 2 && options.colin_sides >= 2;
}

std::vector<long double> distribution(const int dice, const int sides) {
    std::vector<long double> dp(static_cast<std::size_t>(dice * sides + 1), 0.0L);
    dp[0] = 1.0L;

    for (int i = 0; i < dice; ++i) {
        std::vector<long double> next(static_cast<std::size_t>(dice * sides + 1), 0.0L);
        for (std::size_t sum = 0; sum < dp.size(); ++sum) {
            const long double ways = dp[sum];
            if (ways == 0.0L) {
                continue;
            }
            for (int face = 1; face <= sides; ++face) {
                next[sum + static_cast<std::size_t>(face)] += ways;
            }
        }
        dp = std::move(next);
    }

    return dp;
}

long double win_probability(const int pete_dice,
                            const int pete_sides,
                            const int colin_dice,
                            const int colin_sides) {
    const std::vector<long double> pete = distribution(pete_dice, pete_sides);
    const std::vector<long double> colin = distribution(colin_dice, colin_sides);

    long double wins = 0.0L;
    long double total = 0.0L;

    for (std::size_t p = 0; p < pete.size(); ++p) {
        if (pete[p] == 0.0L) {
            continue;
        }
        for (std::size_t c = 0; c < colin.size(); ++c) {
            if (colin[c] == 0.0L) {
                continue;
            }
            const long double ways = pete[p] * colin[c];
            total += ways;
            if (p > c) {
                wins += ways;
            }
        }
    }

    return wins / total;
}

bool run_checkpoints() {
    if (std::abs(win_probability(1, 4, 1, 6) - 0.25L) > 1e-15L) {
        std::cerr << "Checkpoint failed for 1d4 vs 1d6" << '\n';
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

    std::cout << std::fixed << std::setprecision(7)
              << static_cast<double>(win_probability(options.pete_dice,
                                                     options.pete_sides,
                                                     options.colin_dice,
                                                     options.colin_sides))
              << '\n';
    return 0;
}
