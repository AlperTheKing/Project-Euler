#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int defects = 20000;
    int chips = 1000000;
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
        if (parse_int_after_prefix(arg, "--defects=", options.defects) ||
            parse_int_after_prefix(arg, "--chips=", options.chips)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.defects >= 0 && options.chips >= 1;
}

double solve(const int defects, const int chips) {
    if (defects <= 2) {
        return 0.0;
    }

    double answer = 1.0;
    for (int i = 0; i <= defects / 2; ++i) {
        double logp = 0.0;

        for (int j = 0; j < i; ++j) {
            logp += std::log((static_cast<double>(defects - 2 * j) *
                              static_cast<double>(defects - 1 - 2 * j)) /
                             (2.0 * static_cast<double>(chips) * static_cast<double>(j + 1)));
        }
        for (int j = 0; j < defects - i; ++j) {
            logp += std::log(1.0 - static_cast<double>(j) / static_cast<double>(chips));
        }

        answer -= std::exp(logp);
    }

    return answer;
}

bool run_checkpoints() {
    const double sample = solve(3, 7);
    if (std::fabs(sample - 0.0204081633) > 1e-10) {
        std::cerr << "Checkpoint failed for p(3,7)" << '\n';
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

    std::cout << std::fixed << std::setprecision(10) << solve(options.defects, options.chips) << '\n';
    return 0;
}
