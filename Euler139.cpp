#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 100000000;
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

i64 solve(const int limit) {
    i64 count = 0;

    for (int m = 2; 2 * m * (m + 1) < limit; ++m) {
        for (int n = 1; n < m; ++n) {
            if (((m - n) & 1) == 0) {
                continue;
            }
            if (std::gcd(m, n) != 1) {
                continue;
            }

            const int a = m * m - n * n;
            const int b = 2 * m * n;
            const int c = m * m + n * n;
            const int diff = std::abs(a - b);
            if (c % diff != 0) {
                continue;
            }

            const int perimeter = a + b + c;
            count += (limit - 1) / perimeter;
        }
    }

    return count;
}

bool run_checkpoints() {
    if (solve(100) != 9) {
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
