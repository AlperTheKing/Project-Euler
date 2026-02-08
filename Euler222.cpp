#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int pipe_radius = 50;
    int min_radius = 30;
    int max_radius = 50;
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
        if (parse_int_after_prefix(arg, "--pipe-radius=", options.pipe_radius) ||
            parse_int_after_prefix(arg, "--min-radius=", options.min_radius) ||
            parse_int_after_prefix(arg, "--max-radius=", options.max_radius)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.min_radius >= 1 && options.max_radius >= options.min_radius &&
           options.pipe_radius >= options.max_radius;
}

double distance_y(const int pipe_radius, const int a, const int b) {
    const double sum = static_cast<double>(a + b);
    const double side = static_cast<double>(2 * pipe_radius) - sum;
    return std::sqrt(sum * sum - side * side);
}

double solve_optimized(const int pipe_radius, const int min_radius, const int max_radius) {
    const int num_balls = max_radius - min_radius + 1;
    if (num_balls == 1) {
        return static_cast<double>(max_radius * 2);
    }

    const int shift = num_balls - 2;
    std::vector<double> cache(static_cast<std::size_t>(1U << shift) * static_cast<std::size_t>(num_balls), -1.0);

    const auto index_of = [&](const int radius) { return radius - min_radius; };

    std::function<double(unsigned int, int)> dfs = [&](const unsigned int mask, const int last_radius) -> double {
        if (mask == 0U) {
            return distance_y(pipe_radius, last_radius, max_radius) + static_cast<double>(max_radius);
        }

        const std::size_t id = static_cast<std::size_t>(mask) +
                               (static_cast<std::size_t>(index_of(last_radius)) << shift);
        if (cache[id] >= 0.0) {
            return cache[id];
        }

        double best = std::numeric_limits<double>::infinity();
        for (int radius = min_radius; radius <= max_radius; ++radius) {
            const unsigned int bit = 1U << index_of(radius);
            if ((mask & bit) == 0U) {
                continue;
            }
            const double cand = distance_y(pipe_radius, radius, last_radius) + dfs(mask & ~bit, radius);
            best = std::min(best, cand);
        }

        cache[id] = best;
        return best;
    };

    unsigned int mask = (1U << num_balls) - 1U;
    mask &= ~(1U << index_of(max_radius));
    mask &= ~(1U << index_of(max_radius - 1));

    const int first = max_radius - 1;
    return static_cast<double>(first) + dfs(mask, first);
}

double solve_bruteforce(const int pipe_radius, std::vector<int> radii) {
    std::sort(radii.begin(), radii.end());
    double best = std::numeric_limits<double>::infinity();
    do {
        double length = static_cast<double>(radii.front());
        for (std::size_t i = 0; i + 1 < radii.size(); ++i) {
            length += distance_y(pipe_radius, radii[i], radii[i + 1]);
        }
        length += static_cast<double>(radii.back());
        best = std::min(best, length);
    } while (std::next_permutation(radii.begin(), radii.end()));
    return best;
}

bool run_checkpoints() {
    {
        std::vector<int> radii{5, 6, 7, 8};
        const double brute = solve_bruteforce(8, radii);
        const double fast = solve_optimized(8, 5, 8);
        if (std::fabs(brute - fast) > 1e-9) {
            std::cerr << "Checkpoint failed for optimized-vs-bruteforce at (8,5..8)" << '\n';
            return false;
        }
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

    const double best = solve_optimized(options.pipe_radius, options.min_radius, options.max_radius);
    std::cout << static_cast<std::uint64_t>(std::llround(1000.0 * best)) << '\n';
    return 0;
}
