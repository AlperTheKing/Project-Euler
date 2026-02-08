#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int depth = 10;
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
        if (parse_int_after_prefix(arg, "--depth=", options.depth)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.depth >= 0;
}

long double fill_gap(const long double k1,
                     const long double k2,
                     const long double k3,
                     const int depth) {
    if (depth == 0) {
        return 0.0L;
    }

    const long double k4 =
        k1 + k2 + k3 + 2.0L * std::sqrt(k1 * k2 + k2 * k3 + k3 * k1);
    const long double area = std::acos(-1.0L) / (k4 * k4);

    return area +
           fill_gap(k1, k2, k4, depth - 1) +
           fill_gap(k1, k3, k4, depth - 1) +
           fill_gap(k2, k3, k4, depth - 1);
}

long double uncovered_ratio(const int depth) {
    const long double pi = std::acos(-1.0L);
    const long double k = 1.0L + 2.0L / std::sqrt(3.0L);
    const long double r = 1.0L / k;

    long double covered = 3.0L * pi * r * r;
    covered += fill_gap(k, k, k, depth);
    covered += 3.0L * fill_gap(-1.0L, k, k, depth);

    return 1.0L - covered / pi;
}

bool run_checkpoints() {
    if (std::abs(uncovered_ratio(0) - 0.3538290724795825L) > 1e-15L) {
        std::cerr << "Checkpoint failed for depth=0" << '\n';
        return false;
    }
    if (std::abs(uncovered_ratio(1) - 0.19813388055841763L) > 1e-15L) {
        std::cerr << "Checkpoint failed for depth=1" << '\n';
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

    std::cout << std::fixed << std::setprecision(8)
              << static_cast<double>(uncovered_ratio(options.depth)) << '\n';
    return 0;
}
