#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 50000000;
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

int solve(const int limit) {
    std::vector<std::uint8_t> count(static_cast<std::size_t>(limit), 0);

    for (int u = 1; u < limit; ++u) {
        const int vmax = (limit - 1) / u;
        const int upper = std::min(vmax, 3 * u - 1);
        if (upper < 1) {
            continue;
        }

        int v = (-u) % 4;
        if (v < 0) {
            v += 4;
        }
        if (v == 0) {
            v = 4;
        }

        for (; v <= upper; v += 4) {
            const int n = u * v;
            if (count[static_cast<std::size_t>(n)] < 2) {
                ++count[static_cast<std::size_t>(n)];
            }
        }
    }

    int answer = 0;
    for (int n = 1; n < limit; ++n) {
        if (count[static_cast<std::size_t>(n)] == 1) {
            ++answer;
        }
    }
    return answer;
}

bool run_checkpoints() {
    if (solve(100) != 25) {
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
