#include <cmath>
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
    return options.limit >= 12;
}

int solve(const int limit) {
    std::vector<int> solutions(static_cast<std::size_t>(limit + 1), 0);

    for (int a = 1; a <= limit / 3; ++a) {
        for (int b = a; b <= (limit - a) / 2; ++b) {
            const int c2 = a * a + b * b;
            const int c = static_cast<int>(std::sqrt(static_cast<double>(c2)));
            if (c * c != c2) {
                continue;
            }
            const int p = a + b + c;
            if (p <= limit) {
                ++solutions[static_cast<std::size_t>(p)];
            }
        }
    }

    int best_p = 0;
    int best_count = -1;
    for (int p = 0; p <= limit; ++p) {
        if (solutions[static_cast<std::size_t>(p)] > best_count) {
            best_count = solutions[static_cast<std::size_t>(p)];
            best_p = p;
        }
    }

    return best_p;
}

bool run_checkpoints() {
    if (solve(120) != 120) {
        std::cerr << "Checkpoint failed for limit=120" << '\n';
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
