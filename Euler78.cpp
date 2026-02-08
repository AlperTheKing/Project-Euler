#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int mod = 1000000;
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
        if (parse_int_after_prefix(arg, "--mod=", options.mod)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.mod >= 2;
}

int solve(const int mod) {
    std::vector<int> p;
    p.reserve(60000U);
    p.push_back(1);  // p(0)

    for (int n = 1;; ++n) {
        i64 value = 0;

        for (int k = 1;; ++k) {
            const int g1 = k * (3 * k - 1) / 2;
            const int g2 = k * (3 * k + 1) / 2;
            if (g1 > n) {
                break;
            }

            const int sign = (k & 1) ? 1 : -1;
            value += sign * p[static_cast<std::size_t>(n - g1)];
            if (g2 <= n) {
                value += sign * p[static_cast<std::size_t>(n - g2)];
            }
        }

        int modded = static_cast<int>(value % mod);
        if (modded < 0) {
            modded += mod;
        }
        p.push_back(modded);

        if (modded == 0) {
            return n;
        }
    }
}

bool run_checkpoints() {
    if (solve(5) != 4) {
        std::cerr << "Checkpoint failed for mod=5" << '\n';
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

    std::cout << solve(options.mod) << '\n';
    return 0;
}
