#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 1000000;
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
    return options.limit >= 1;
}

std::vector<int> compute_grundy(const int n) {
    std::vector<int> grundy(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> seen(2048, -1);

    for (int len = 1; len <= n; ++len) {
        for (int left = 0; left <= len - 2 - left; ++left) {
            const int value = grundy[static_cast<std::size_t>(left)] ^
                              grundy[static_cast<std::size_t>(len - 2 - left)];
            seen[static_cast<std::size_t>(value)] = len;
        }
        while (seen[static_cast<std::size_t>(grundy[static_cast<std::size_t>(len)])] == len) {
            ++grundy[static_cast<std::size_t>(len)];
        }
    }
    return grundy;
}

int solve(const int limit) {
    const int prefix = 100;
    const int period = 34;
    const int max_precompute = std::max(prefix + period, 1200);
    const std::vector<int> g = compute_grundy(max_precompute);

    int wins = 0;
    for (int len = 1; len <= limit; ++len) {
        const int sg = (len <= prefix) ? g[static_cast<std::size_t>(len)]
                                       : g[static_cast<std::size_t>((len - prefix) % period + prefix)];
        wins += (sg != 0);
    }
    return wins;
}

bool run_checkpoints() {
    const std::vector<int> g = compute_grundy(1200);
    for (int i = 100; i <= 1200; ++i) {
        if (g[static_cast<std::size_t>(i)] != g[static_cast<std::size_t>((i - 100) % 34 + 100)]) {
            std::cerr << "Checkpoint failed for Grundy periodicity" << '\n';
            return false;
        }
    }
    if (solve(5) != 3) {
        std::cerr << "Checkpoint failed for [1,5]" << '\n';
        return false;
    }
    if (solve(50) != 40) {
        std::cerr << "Checkpoint failed for [1,50]" << '\n';
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
