#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int shots = 50;
    int target_score = 20;
    long double target_probability = 0.02L;
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

bool parse_ld_after_prefix(const std::string& arg, const std::string& prefix, long double& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    char* end_ptr = nullptr;
    value = std::strtold(tail.c_str(), &end_ptr);
    return end_ptr != nullptr && *end_ptr == '\0';
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--shots=", options.shots) ||
            parse_int_after_prefix(arg, "--target-score=", options.target_score) ||
            parse_ld_after_prefix(arg, "--target-probability=", options.target_probability)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.shots >= 1 && options.target_score >= 0 &&
           options.target_score <= options.shots &&
           options.target_probability > 0.0L && options.target_probability < 1.0L;
}

long double exact_score_probability(const int shots, const int target_score, const long double q) {
    std::array<long double, 64> dp{};
    dp.fill(0.0L);
    dp[0] = 1.0L;
    for (int x = 1; x <= shots; ++x) {
        const long double hit = 1.0L - static_cast<long double>(x) / q;
        const long double miss = 1.0L - hit;
        for (int s = x; s >= 1; --s) {
            dp[static_cast<std::size_t>(s)] =
                dp[static_cast<std::size_t>(s)] * miss +
                dp[static_cast<std::size_t>(s - 1)] * hit;
        }
        dp[0] *= miss;
    }
    return dp[static_cast<std::size_t>(target_score)];
}

long double solve(const Options& options) {
    long double low = 50.0L;
    long double high = 100.0L;
    while (exact_score_probability(options.shots, options.target_score, high) > options.target_probability) {
        high *= 2.0L;
    }

    for (int it = 0; it < 220; ++it) {
        const long double mid = (low + high) / 2.0L;
        if (exact_score_probability(options.shots, options.target_score, mid) > options.target_probability) {
            low = mid;
        } else {
            high = mid;
        }
    }
    return (low + high) / 2.0L;
}

bool run_checkpoints() {
    const long double p52 = exact_score_probability(50, 20, 52.0L);
    const long double p53 = exact_score_probability(50, 20, 53.0L);
    if (!(p52 > 0.02L && p53 < 0.02L)) {
        std::cerr << "Checkpoint failed for bracketing around q=52..53" << '\n';
        return false;
    }
    const Options options;
    const long double q = solve(options);
    const long double prob = exact_score_probability(options.shots, options.target_score, q);
    if (std::fabsl(prob - options.target_probability) > 1e-13L) {
        std::cerr << "Checkpoint failed for solved q probability consistency" << '\n';
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
    const long double answer = solve(options);
    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(answer) << '\n';
    return 0;
}
