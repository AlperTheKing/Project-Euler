#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int from = 2;
    int to = 15;
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
        if (parse_int_after_prefix(arg, "--from=", options.from) ||
            parse_int_after_prefix(arg, "--to=", options.to)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.from >= 2 && options.to >= options.from;
}

i64 P_floor(const int m) {
    long double log_sum = 0.0L;
    for (int i = 1; i <= m; ++i) {
        const long double xi = 2.0L * static_cast<long double>(i) / static_cast<long double>(m + 1);
        log_sum += static_cast<long double>(i) * std::log(xi);
    }

    const long double value = std::exp(log_sum);
    return static_cast<i64>(std::floor(value + 1e-12L));
}

i64 solve(const int from, const int to) {
    i64 sum = 0;
    for (int m = from; m <= to; ++m) {
        sum += P_floor(m);
    }
    return sum;
}

bool run_checkpoints() {
    if (P_floor(10) != 4112LL) {
        std::cerr << "Checkpoint failed for P(10)" << '\n';
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

    std::cout << solve(options.from, options.to) << '\n';
    return 0;
}
