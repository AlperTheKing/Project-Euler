#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int grid_size = 20;
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

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--grid-size=", options.grid_size)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.grid_size >= 0;
}

u64 binomial(const int n, int k) {
    k = std::min(k, n - k);
    u128 value = 1;

    for (int i = 1; i <= k; ++i) {
        value = (value * static_cast<u128>(n - k + i)) / static_cast<u128>(i);
    }

    if (value > static_cast<u128>(std::numeric_limits<u64>::max())) {
        throw std::overflow_error("Result overflow");
    }
    return static_cast<u64>(value);
}

u64 solve(const int grid_size) {
    return binomial(2 * grid_size, grid_size);
}

bool run_checkpoints() {
    if (solve(2) != 6ULL) {
        std::cerr << "Checkpoint failed for grid_size=2" << '\n';
        return false;
    }
    if (solve(1) != 2ULL) {
        std::cerr << "Checkpoint failed for grid_size=1" << '\n';
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

    std::cout << solve(options.grid_size) << '\n';
    return 0;
}
