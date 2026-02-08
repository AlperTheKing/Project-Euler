#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int limit = 64000000;
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

bool is_perfect_square(const u64 x) {
    const u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    return r * r == x || (r + 1) * (r + 1) == x;
}

u64 solve(const int limit) {
    std::vector<u64> sigma2(static_cast<std::size_t>(limit + 1), 1ULL);
    sigma2[0] = 0;

    for (int d = 2; d <= limit; ++d) {
        const u64 d2 = static_cast<u64>(d) * static_cast<u64>(d);
        for (int m = d; m <= limit; m += d) {
            sigma2[static_cast<std::size_t>(m)] += d2;
        }
    }

    u64 sum = 0;
    for (int n = 1; n <= limit; ++n) {
        if (is_perfect_square(sigma2[static_cast<std::size_t>(n)])) {
            sum += static_cast<u64>(n);
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(10) != 1ULL) {
        std::cerr << "Checkpoint failed for limit=10" << '\n';
        return false;
    }
    if (solve(1000) != 1304ULL) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
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
