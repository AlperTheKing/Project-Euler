#include <cmath>
#include <cstdint>
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

double area_A(const int p, const int q) {
    const double dp = static_cast<double>(p);
    const double delta = std::acos(-1.0) / dp;
    const double sd = std::sin(delta);
    const double cd = std::cos(delta);
    const double tan_delta = sd / cd;

    double sum_diff = 0.0;
    if (q >= 2) {
        int k = 2;
        double sin_k = sd;
        double cos_k = cd;
        double sign = 1.0;
        double corr = 0.0;  // Kahan compensation.
        while (k <= q) {
            const double sin_next = sin_k * cd + cos_k * sd;
            const double cos_next = cos_k * cd - sin_k * sd;
            const double diff = sd / (cos_k * cos_next);
            const double term = sign * diff;

            const double y = term - corr;
            const double t = sum_diff + y;
            corr = (t - sum_diff) - y;
            sum_diff = t;

            sin_k = sin_next;
            cos_k = cos_next;
            sign = -sign;
            ++k;
        }
    }

    const double inner = sum_diff - tan_delta;
    const double outside_sign = ((q % 2) == 0) ? 1.0 : -1.0;
    return dp * outside_sign * inner;
}

std::vector<std::uint64_t> fibonacci_sequence(const int count) {
    std::vector<std::uint64_t> fibs;
    fibs.reserve(static_cast<std::size_t>(count));
    fibs.push_back(1ULL);
    fibs.push_back(1ULL);
    for (int i = 2; i < count; ++i) {
        fibs.push_back(fibs[static_cast<std::size_t>(i - 1)] + fibs[static_cast<std::size_t>(i - 2)]);
    }
    return fibs;
}

double solve() {
    const std::vector<std::uint64_t> fibs = fibonacci_sequence(36);  // F1..F36
    double sum = 0.0;
    double corr = 0.0;
    for (int n = 3; n <= 34; ++n) {
        const int p = static_cast<int>(fibs[static_cast<std::size_t>(n)]);
        const int q = static_cast<int>(fibs[static_cast<std::size_t>(n - 2)]);
        const double term = area_A(p, q);
        const double y = term - corr;
        const double t = sum + y;
        corr = (t - sum) - y;
        sum = t;
    }
    return sum;
}

bool run_checkpoints() {
    const double a1 = area_A(8, 3);
    if (std::fabs(a1 - 9.9411254970) > 5e-11) {
        std::cerr << "Checkpoint failed: A(8,3)" << '\n';
        return false;
    }
    const double a2 = area_A(130021, 50008);
    if (std::fabs(a2 - 10.9210371479) > 5e-11) {
        std::cerr << "Checkpoint failed: A(130021,50008)" << '\n';
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
    std::cout << std::fixed << std::setprecision(10) << solve() << '\n';
    return 0;
}
