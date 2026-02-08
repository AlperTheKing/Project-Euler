#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    long double eps = 1e-12L;
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
    return value > 0;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_ld_after_prefix(arg, "--eps=", options.eps)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.eps > 0;
}

long double s_fn(long double x) {
    x = std::fabsl(x);
    const long double nearest = std::floor(x + 0.5L);
    return std::fabsl(x - nearest);
}

long double blancmange(long double x) {
    long double y = 0.0L;
    long double scale = 1.0L;

    for (int i = 0; i < 60; ++i) {
        y += s_fn(scale * x) / scale;
        scale *= 2.0L;
    }

    return y;
}

long double overlap_height(const long double x) {
    const long double dx = x - 0.25L;
    const long double inside = 0.0625L - dx * dx;
    if (inside <= 0.0L) {
        return 0.0L;
    }

    const long double dy = std::sqrt(inside);
    const long double circle_low = 0.5L - dy;
    const long double circle_high = 0.5L + dy;

    const long double under_curve = blancmange(x);
    const long double lo = std::max(0.0L, circle_low);
    const long double hi = std::min(circle_high, under_curve);
    if (hi <= lo) {
        return 0.0L;
    }
    return hi - lo;
}

long double solve(const long double eps) {
    const long double l = 0.0L;
    const long double r = 0.5L;
    int steps = 1 << 18;
    while (steps < (1 << 22) && (1.0L / static_cast<long double>(steps * steps)) > eps) {
        steps <<= 1;
    }

    const long double h = (r - l) / static_cast<long double>(steps);
    long double odd_sum = 0.0L;
    long double even_sum = 0.0L;

    for (int i = 1; i < steps; ++i) {
        const long double x = l + h * static_cast<long double>(i);
        const long double fx = overlap_height(x);
        if (i & 1) {
            odd_sum += fx;
        } else {
            even_sum += fx;
        }
    }

    const long double f0 = overlap_height(l);
    const long double fn = overlap_height(r);
    return h * (f0 + fn + 4.0L * odd_sum + 2.0L * even_sum) / 3.0L;
}

bool run_checkpoints() {
    if (std::fabsl(blancmange(0.0L) - 0.0L) > 1e-15L) {
        std::cerr << "Checkpoint failed for blancmange(0)" << '\n';
        return false;
    }
    if (std::fabsl(blancmange(0.5L) - 0.5L) > 1e-12L) {
        std::cerr << "Checkpoint failed for blancmange(0.5)" << '\n';
        return false;
    }

    const long double area = solve(1e-12L);
    if (area <= 0.0L || area >= 0.5L) {
        std::cerr << "Checkpoint failed for area range" << '\n';
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

    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(solve(options.eps)) << '\n';
    return 0;
}
