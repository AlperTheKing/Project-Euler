#include <algorithm>
#include <iostream>
#include <string>

namespace {

struct Options {
    int max_multiplier = 6;
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
        if (parse_int_after_prefix(arg, "--max-multiplier=", options.max_multiplier)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.max_multiplier >= 2;
}

std::string sorted_digits(int value) {
    std::string s = std::to_string(value);
    std::sort(s.begin(), s.end());
    return s;
}

int solve(const int max_multiplier) {
    for (int x = 1;; ++x) {
        const std::string signature = sorted_digits(x);
        bool ok = true;
        for (int k = 2; k <= max_multiplier; ++k) {
            if (sorted_digits(k * x) != signature) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return x;
        }
    }
}

bool run_checkpoints() {
    if (sorted_digits(125874) != sorted_digits(251748)) {
        std::cerr << "Checkpoint failed for 125874" << '\n';
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

    std::cout << solve(options.max_multiplier) << '\n';
    return 0;
}
