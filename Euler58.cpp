#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    int threshold_percent = 10;
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
        if (parse_int_after_prefix(arg, "--threshold-percent=", options.threshold_percent)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.threshold_percent >= 1 && options.threshold_percent <= 100;
}

bool is_prime(const i64 n) {
    if (n < 2) {
        return false;
    }
    if ((n % 2) == 0) {
        return n == 2;
    }
    for (i64 p = 3; p <= n / p; p += 2) {
        if ((n % p) == 0) {
            return false;
        }
    }
    return true;
}

int solve(const int threshold_percent) {
    i64 prime_count = 0;
    i64 diagonal_count = 1;

    for (int layer = 1;; ++layer) {
        const i64 side = 2LL * layer + 1LL;
        const i64 square = side * side;
        const i64 step = side - 1;

        if (is_prime(square - step)) {
            ++prime_count;
        }
        if (is_prime(square - 2 * step)) {
            ++prime_count;
        }
        if (is_prime(square - 3 * step)) {
            ++prime_count;
        }

        diagonal_count += 4;

        if (prime_count * 100 < diagonal_count * threshold_percent) {
            return static_cast<int>(side);
        }
    }
}

bool run_checkpoints() {
    if (solve(62) != 3) {
        std::cerr << "Checkpoint failed for threshold 62%" << '\n';
        return false;
    }
    if (solve(60) != 5) {
        std::cerr << "Checkpoint failed for threshold 60%" << '\n';
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

    std::cout << solve(options.threshold_percent) << '\n';
    return 0;
}
