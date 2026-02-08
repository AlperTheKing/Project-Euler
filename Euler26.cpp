#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 1000;
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
    for (const char c : tail) {
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

    return options.limit > 2;
}

int recurring_cycle_length(const int denominator) {
    std::vector<int> first_seen(static_cast<std::size_t>(denominator), -1);

    int remainder = 1 % denominator;
    int position = 0;

    while (remainder != 0 && first_seen[static_cast<std::size_t>(remainder)] == -1) {
        first_seen[static_cast<std::size_t>(remainder)] = position;
        remainder = (remainder * 10) % denominator;
        ++position;
    }

    if (remainder == 0) {
        return 0;
    }
    return position - first_seen[static_cast<std::size_t>(remainder)];
}

int solve(const int limit) {
    int best_d = 2;
    int best_len = 0;

    for (int d = 2; d < limit; ++d) {
        const int cycle = recurring_cycle_length(d);
        if (cycle > best_len) {
            best_len = cycle;
            best_d = d;
        }
    }

    return best_d;
}

bool run_checkpoints() {
    if (recurring_cycle_length(7) != 6) {
        std::cerr << "Checkpoint failed for denominator=7" << '\n';
        return false;
    }
    if (solve(10) != 7) {
        std::cerr << "Checkpoint failed for limit=10" << '\n';
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
