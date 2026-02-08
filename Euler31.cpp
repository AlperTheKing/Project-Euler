#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int target = 200;
    bool run_checkpoints = true;
};

constexpr std::array<int, 8> kCoins = {1, 2, 5, 10, 20, 50, 100, 200};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
        if (parsed > static_cast<i64>(std::numeric_limits<int>::max())) {
            return false;
        }
    }

    value = static_cast<int>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--target=", options.target)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.target >= 0;
}

i64 solve(const int target) {
    std::vector<i64> ways(static_cast<std::size_t>(target + 1), 0LL);
    ways[0] = 1;

    for (const int coin : kCoins) {
        for (int amount = coin; amount <= target; ++amount) {
            ways[static_cast<std::size_t>(amount)] += ways[static_cast<std::size_t>(amount - coin)];
        }
    }

    return ways[static_cast<std::size_t>(target)];
}

bool run_checkpoints() {
    if (solve(5) != 4LL) {
        std::cerr << "Checkpoint failed for target=5" << '\n';
        return false;
    }
    if (solve(10) != 11LL) {
        std::cerr << "Checkpoint failed for target=10" << '\n';
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

    std::cout << solve(options.target) << '\n';
    return 0;
}
