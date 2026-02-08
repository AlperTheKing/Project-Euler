#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int iterations = 100000;
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
        if (parse_int_after_prefix(arg, "--iterations=", options.iterations)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.iterations >= 1;
}

double f(const double x) {
    return std::floor(std::pow(2.0, 30.403243784 - x * x)) * 1e-9;
}

double solve(const int iterations) {
    double u = -1.0;
    for (int i = 0; i < iterations; ++i) {
        u = f(u);
    }
    const double v = f(u);
    return u + v;
}

bool run_checkpoints() {
    if (std::abs(f(-1.0) - 0.7100000000) > 1e-12) {
        std::cerr << "Checkpoint failed for first sequence value" << '\n';
        return false;
    }
    if (std::abs(solve(10000) - 1.710637717) > 1e-12) {
        std::cerr << "Checkpoint failed for converged two-cycle sum" << '\n';
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

    std::cout << std::fixed << std::setprecision(9) << solve(options.iterations) << '\n';
    return 0;
}
