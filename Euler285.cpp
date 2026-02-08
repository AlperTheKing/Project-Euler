#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int max_k = 100000;
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
        if (parse_int_after_prefix(arg, "--max-k=", options.max_k)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_k >= 1;
}

long double area_inside_square_and_circle(const long double r) {
    static const long double kSqrt2 = std::sqrt(2.0L);
    if (r <= kSqrt2) {
        return 0.0L;
    }
    const long double t = std::sqrt(r * r - 1.0L);
    const long double pi = std::acos(-1.0L);
    return r * r * (pi / 4.0L - std::asin(1.0L / r)) - t + 1.0L;
}

long double solve(const int max_k) {
    long double expected = 0.0L;
    for (int k = 1; k <= max_k; ++k) {
        const long double outer = area_inside_square_and_circle(static_cast<long double>(k) + 0.5L);
        const long double inner = area_inside_square_and_circle(static_cast<long double>(k) - 0.5L);
        expected += (outer - inner) / static_cast<long double>(k);
    }
    return expected;
}

bool run_checkpoints() {
    const long double sample = solve(10);
    const long double rounded_sample = std::round(sample * 100000.0L) / 100000.0L;
    if (std::fabsl(rounded_sample - 10.20914L) > 1e-12L) {
        std::cerr << "Checkpoint failed for 10-turn sample" << '\n';
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
    const long double answer = solve(options.max_k);
    std::cout << std::fixed << std::setprecision(5) << static_cast<double>(answer) << '\n';
    return 0;
}
