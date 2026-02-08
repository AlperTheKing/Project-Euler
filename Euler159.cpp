#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

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

    return options.limit >= 2;
}

int digital_root(const int n) {
    return (n == 0) ? 0 : 1 + (n - 1) % 9;
}

std::vector<int> compute_mdrs(const int limit) {
    std::vector<int> mdrs(static_cast<std::size_t>(limit), 0);
    for (int n = 1; n < limit; ++n) {
        mdrs[static_cast<std::size_t>(n)] = digital_root(n);
    }

    for (int a = 2; a < limit; ++a) {
        for (int b = 2, p = a * 2; p < limit; ++b, p += a) {
            const int cand = mdrs[static_cast<std::size_t>(a)] + mdrs[static_cast<std::size_t>(b)];
            if (cand > mdrs[static_cast<std::size_t>(p)]) {
                mdrs[static_cast<std::size_t>(p)] = cand;
            }
        }
    }

    return mdrs;
}

i64 solve(const int limit) {
    const std::vector<int> mdrs = compute_mdrs(limit);

    i64 sum = 0;
    for (int n = 2; n < limit; ++n) {
        sum += mdrs[static_cast<std::size_t>(n)];
    }
    return sum;
}

bool run_checkpoints() {
    const std::vector<int> mdrs = compute_mdrs(100);
    if (mdrs[24] != 11) {
        std::cerr << "Checkpoint failed for mdrs(24)" << '\n';
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
