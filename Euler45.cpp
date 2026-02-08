#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

bool is_pentagonal(const i64 x) {
    if (x <= 0) {
        return false;
    }
    const long double disc = 1.0L + 24.0L * static_cast<long double>(x);
    const long double root = std::sqrt(disc);
    const long double n = (1.0L + root) / 6.0L;
    const i64 k = static_cast<i64>(std::llround(n));
    return k > 0 && (k * (3 * k - 1) / 2 == x);
}

i64 solve() {
    for (i64 n = 144;; ++n) {
        const i64 hexagonal = n * (2 * n - 1);
        if (is_pentagonal(hexagonal)) {
            return hexagonal;
        }
    }
}

bool run_checkpoints() {
    if (!is_pentagonal(40755)) {
        std::cerr << "Checkpoint failed for 40755" << '\n';
        return false;
    }
    if (is_pentagonal(40756)) {
        std::cerr << "Checkpoint failed for 40756" << '\n';
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

    std::cout << solve() << '\n';
    return 0;
}
