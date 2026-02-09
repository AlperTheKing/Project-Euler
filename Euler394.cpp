#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    long double x = 40.0L;
    bool run_checkpoints = true;
};

bool parse_ld_after_prefix(const std::string& arg, const std::string& prefix, long double& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = std::stold(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_ld_after_prefix(arg, "--x=", options.x)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.x >= 1.0L;
}

long double expected_steps(const long double x) {
    // Closed form from the renewal equation:
    // E(x) = 1 + 2 * integral_0^1 (1-w) E(xw) dw (with E(u)=0 for u<1),
    // yielding E(x) = 7/9 + (2/3) ln x + 2/(9 x^3).
    return 7.0L / 9.0L + (2.0L / 3.0L) * std::log(x) + 2.0L / (9.0L * x * x * x);
}

bool nearly_equal(const long double a, const long double b, const long double rel_tol = 1e-12L) {
    const long double scale = std::max(std::fabsl(a), std::fabsl(b));
    if (scale == 0.0L) {
        return true;
    }
    return std::fabsl(a - b) <= rel_tol * scale;
}

bool run_checkpoints() {
    if (!nearly_equal(expected_steps(1.0L), 1.0L, 1e-14L)) {
        std::cerr << "Checkpoint failed: E(1)\n";
        return false;
    }

    if (!nearly_equal(expected_steps(2.0L), 1.2676536759L, 1e-10L)) {
        std::cerr << "Checkpoint failed: E(2)\n";
        return false;
    }

    if (!nearly_equal(expected_steps(7.5L), 2.1215732071L, 1e-10L)) {
        std::cerr << "Checkpoint failed: E(7.5)\n";
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

    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(expected_steps(options.x)) << '\n';
    return 0;
}
