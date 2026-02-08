#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    i64 min_total = 1000000000000LL;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--min-total=", options.min_total)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.min_total >= 1;
}

bool valid_probability(const i64 blue, const i64 total) {
    return 2LL * blue * (blue - 1LL) == total * (total - 1LL);
}

i64 solve(const i64 min_total) {
    i64 blue = 15;
    i64 total = 21;

    while (total <= min_total) {
        const i64 next_blue = 3LL * blue + 2LL * total - 2LL;
        const i64 next_total = 4LL * blue + 3LL * total - 3LL;
        blue = next_blue;
        total = next_total;
    }

    return blue;
}

bool run_checkpoints() {
    if (!valid_probability(15, 21)) {
        std::cerr << "Checkpoint failed for (15,21) probability" << '\n';
        return false;
    }
    if (solve(21) != 85) {
        std::cerr << "Checkpoint failed for limit 21" << '\n';
        return false;
    }
    if (solve(120) != 493) {
        std::cerr << "Checkpoint failed for limit 120" << '\n';
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

    std::cout << solve(options.min_total) << '\n';
    return 0;
}
