#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int n = 5000;
    long double target = -600000000000.0L;
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
        if (parse_int_after_prefix(arg, "--n=", options.n) ||
            parse_ld_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

long double series_sum(const int n, const long double r) {
    long double sum = 0.0L;
    long double power = 1.0L;
    for (int k = 1; k <= n; ++k) {
        sum += (900.0L - 3.0L * static_cast<long double>(k)) * power;
        power *= r;
    }
    return sum;
}

long double solve(const int n, const long double target) {
    long double lo = 1.0L;
    long double hi = 1.2L;

    for (int it = 0; it < 200; ++it) {
        const long double mid = (lo + hi) / 2.0L;
        const long double s = series_sum(n, mid);
        if (s > target) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return (lo + hi) / 2.0L;
}

bool run_checkpoints() {
    if (std::abs(series_sum(5000, 1.0L) + 33007500.0L) > 1e-6L) {
        std::cerr << "Checkpoint failed for s(5000) at r=1" << '\n';
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

    std::cout << std::fixed << std::setprecision(12)
              << static_cast<double>(solve(options.n, options.target)) << '\n';
    return 0;
}
