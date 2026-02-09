#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int colors = 7;
    int balls_per_color = 10;
    int draws = 20;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--colors=", options.colors) ||
            parse_int_after_prefix(arg, "--balls-per-color=", options.balls_per_color) ||
            parse_int_after_prefix(arg, "--draws=", options.draws)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.colors >= 1 && options.balls_per_color >= 1 && options.draws >= 0 &&
           options.draws <= options.colors * options.balls_per_color;
}

long double choose_ratio(const int a, const int b, const int draws) {
    // Computes C(a,draws)/C(b,draws) as a stable product.
    if (draws == 0) {
        return 1.0L;
    }
    long double ratio = 1.0L;
    for (int i = 0; i < draws; ++i) {
        ratio *= static_cast<long double>(a - i) / static_cast<long double>(b - i);
    }
    return ratio;
}

long double expected_distinct_colors(const int colors, const int balls_per_color, const int draws) {
    const int total = colors * balls_per_color;
    const int without_one_color = (colors - 1) * balls_per_color;
    const long double miss_one = choose_ratio(without_one_color, total, draws);
    return static_cast<long double>(colors) * (1.0L - miss_one);
}

long double brute_tiny_expectation(const int colors, const int balls_per_color, const int draws) {
    const int total = colors * balls_per_color;
    std::vector<int> bag;
    bag.reserve(static_cast<std::size_t>(total));
    for (int c = 0; c < colors; ++c) {
        for (int t = 0; t < balls_per_color; ++t) {
            bag.push_back(c);
        }
    }

    std::vector<int> pick(draws);
    long double sum = 0.0L;
    std::uint64_t ways = 0ULL;

    std::function<void(int, int)> dfs = [&](int idx, int start) {
        if (idx == draws) {
            std::vector<unsigned char> seen(static_cast<std::size_t>(colors), 0U);
            int distinct = 0;
            for (int x : pick) {
                if (seen[static_cast<std::size_t>(x)] == 0U) {
                    seen[static_cast<std::size_t>(x)] = 1U;
                    ++distinct;
                }
            }
            sum += static_cast<long double>(distinct);
            ++ways;
            return;
        }
        for (int i = start; i <= total - (draws - idx); ++i) {
            pick[static_cast<std::size_t>(idx)] = bag[static_cast<std::size_t>(i)];
            dfs(idx + 1, i + 1);
        }
    };
    dfs(0, 0);
    return sum / static_cast<long double>(ways);
}

bool run_checkpoints() {
    const long double exact = expected_distinct_colors(3, 2, 2);
    const long double brute = brute_tiny_expectation(3, 2, 2);
    if (std::fabsl(exact - brute) > 1e-15L) {
        std::cerr << "Checkpoint failed: tiny hypergeometric expectation cross-check" << '\n';
        return false;
    }
    if (std::fabsl(expected_distinct_colors(1, 10, 0) - 0.0L) > 1e-18L) {
        std::cerr << "Checkpoint failed: zero-draw edge case" << '\n';
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
    const long double ans = expected_distinct_colors(options.colors, options.balls_per_color, options.draws);
    std::cout << std::fixed << std::setprecision(9) << static_cast<double>(ans) << '\n';
    return 0;
}
