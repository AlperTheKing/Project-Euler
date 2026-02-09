#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

struct Stats {
    long double mean{0.0L};
    long double variance{0.0L};
};

Stats uniform_die_stats(const int sides) {
    const long double s = static_cast<long double>(sides);
    return Stats{(s + 1.0L) / 2.0L, (s * s - 1.0L) / 12.0L};
}

Stats transition_sum_of_random_dice(const Stats& count_stats, const int sides) {
    const Stats per_die = uniform_die_stats(sides);
    const long double mean =
        per_die.mean * count_stats.mean;
    const long double variance =
        per_die.variance * count_stats.mean +
        per_die.mean * per_die.mean * count_stats.variance;
    return Stats{mean, variance};
}

Stats brute_one_layer_uniform_count(const int max_count, const int sides) {
    // N is uniform on {1..max_count}; then S is sum of N fair `sides`-sided dice.
    const long double p_count = 1.0L / static_cast<long double>(max_count);
    const int max_sum = max_count * sides;
    std::vector<long double> total_dist(static_cast<std::size_t>(max_sum + 1), 0.0L);

    for (int n = 1; n <= max_count; ++n) {
        std::vector<long double> dist(1, 1.0L);  // sum=0 before rolling.
        for (int roll = 0; roll < n; ++roll) {
            std::vector<long double> next(static_cast<std::size_t>(dist.size() + sides), 0.0L);
            for (std::size_t sum = 0; sum < dist.size(); ++sum) {
                const long double base = dist[sum] / static_cast<long double>(sides);
                for (int face = 1; face <= sides; ++face) {
                    next[sum + static_cast<std::size_t>(face)] += base;
                }
            }
            dist.swap(next);
        }

        for (std::size_t sum = 0; sum < dist.size(); ++sum) {
            total_dist[sum] += p_count * dist[sum];
        }
    }

    long double mean = 0.0L;
    for (std::size_t sum = 0; sum < total_dist.size(); ++sum) {
        mean += total_dist[sum] * static_cast<long double>(sum);
    }

    long double variance = 0.0L;
    for (std::size_t sum = 0; sum < total_dist.size(); ++sum) {
        const long double diff = static_cast<long double>(sum) - mean;
        variance += total_dist[sum] * diff * diff;
    }

    return Stats{mean, variance};
}

bool nearly_equal(const long double a, const long double b, const long double rel_tol = 1e-14L) {
    const long double scale = std::max(std::fabsl(a), std::fabsl(b));
    if (scale == 0.0L) {
        return true;
    }
    return std::fabsl(a - b) <= rel_tol * scale;
}

bool run_checkpoints() {
    const Stats t = uniform_die_stats(4);
    if (!nearly_equal(t.mean, 2.5L) || !nearly_equal(t.variance, 1.25L)) {
        std::cerr << "Checkpoint failed for T\n";
        return false;
    }

    const Stats formula_c = transition_sum_of_random_dice(t, 6);
    const Stats brute_c = brute_one_layer_uniform_count(4, 6);
    if (!nearly_equal(formula_c.mean, brute_c.mean) ||
        !nearly_equal(formula_c.variance, brute_c.variance)) {
        std::cerr << "Checkpoint failed for C (formula vs brute distribution)\n";
        return false;
    }

    return true;
}

long double solve_variance() {
    const Stats t = uniform_die_stats(4);
    const Stats c = transition_sum_of_random_dice(t, 6);
    const Stats o = transition_sum_of_random_dice(c, 8);
    const Stats d = transition_sum_of_random_dice(o, 12);
    const Stats i = transition_sum_of_random_dice(d, 20);
    return i.variance;
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

    std::cout << std::fixed << std::setprecision(4) << static_cast<double>(solve_variance()) << '\n';
    return 0;
}
