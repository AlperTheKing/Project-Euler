#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 10000000;
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

    return options.limit >= 3;
}

int solve(const int limit) {
    std::vector<std::uint16_t> div_count(static_cast<std::size_t>(limit + 2), 0);

    for (int d = 1; d <= limit + 1; ++d) {
        for (int m = d; m <= limit + 1; m += d) {
            ++div_count[static_cast<std::size_t>(m)];
        }
    }

    int answer = 0;
    for (int n = 2; n < limit; ++n) {
        if (div_count[static_cast<std::size_t>(n)] == div_count[static_cast<std::size_t>(n + 1)]) {
            ++answer;
        }
    }

    return answer;
}

bool run_checkpoints() {
    if (solve(100) != 15) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
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
