#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int tosses = 1000;
    long double target = 1.0e9L;
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
        if (parse_int_after_prefix(arg, "--tosses=", options.tosses) ||
            parse_ld_after_prefix(arg, "--target=", options.target)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.tosses >= 1 && options.target > 1.0L;
}

long double binomial_tail_probability(const int n, const int k_min) {
    if (k_min <= 0) {
        return 1.0L;
    }
    if (k_min > n) {
        return 0.0L;
    }

    std::vector<long double> prob(static_cast<std::size_t>(n + 1), 0.0L);
    prob[0] = std::ldexp(1.0L, -n);  // 2^-n
    for (int k = 0; k < n; ++k) {
        prob[static_cast<std::size_t>(k + 1)] =
            prob[static_cast<std::size_t>(k)] * static_cast<long double>(n - k) /
            static_cast<long double>(k + 1);
    }

    long double sum = 0.0L;
    for (int k = k_min; k <= n; ++k) {
        sum += prob[static_cast<std::size_t>(k)];
    }
    return sum;
}

int minimum_wins_needed(const int n, const long double target) {
    const long double ln_target = std::log(target);
    for (int k = 0; k <= n; ++k) {
        long double best_gain = -1.0e300L;
        if (k == n) {
            best_gain = static_cast<long double>(n) * std::log(3.0L);
        } else {
            const long double f_star =
                (3.0L * static_cast<long double>(k) - static_cast<long double>(n)) /
                (2.0L * static_cast<long double>(n));
            if (f_star > 0.0L && f_star < 1.0L) {
                best_gain =
                    static_cast<long double>(k) * std::log(1.0L + 2.0L * f_star) +
                    static_cast<long double>(n - k) * std::log(1.0L - f_star);
            } else if (f_star <= 0.0L) {
                best_gain = 0.0L;
            }
        }
        if (best_gain >= ln_target) {
            return k;
        }
    }
    return n + 1;
}

long double solve_probability(const int tosses, const long double target) {
    const int k_min = minimum_wins_needed(tosses, target);
    return binomial_tail_probability(tosses, k_min);
}

bool run_checkpoints() {
    // For two tosses and target 2.0, requiring exactly two wins gives probability 1/4.
    {
        const long double p = solve_probability(2, 2.0L);
        if (std::fabsl(p - 0.25L) > 1e-15L) {
            std::cerr << "Checkpoint failed for tosses=2 target=2.0" << '\n';
            return false;
        }
    }
    // For four tosses and target 4.0, at least three wins are required at optimum.
    {
        const long double p = solve_probability(4, 4.0L);
        const long double expected = 5.0L / 16.0L;
        if (std::fabsl(p - expected) > 1e-15L) {
            std::cerr << "Checkpoint failed for tosses=4 target=4.0" << '\n';
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
    const long double ans = solve_probability(options.tosses, options.target);
    std::cout << std::fixed << std::setprecision(12) << static_cast<double>(ans) << '\n';
    return 0;
}
