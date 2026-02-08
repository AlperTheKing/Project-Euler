#include <cmath>
#include <iostream>
#include <string>

namespace {

struct Options {
    int limit = 10000;
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
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 2;
}

int period_length(const int n) {
    const int a0 = static_cast<int>(std::sqrt(n));
    if (a0 * a0 == n) {
        return 0;
    }

    int m = 0;
    int d = 1;
    int a = a0;
    int period = 0;

    do {
        m = d * a - m;
        d = (n - m * m) / d;
        a = (a0 + m) / d;
        ++period;
    } while (a != 2 * a0);

    return period;
}

int solve(const int limit) {
    int count = 0;
    for (int n = 2; n <= limit; ++n) {
        if ((period_length(n) & 1) == 1) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(13) != 4) {
        std::cerr << "Checkpoint failed for limit=13" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
