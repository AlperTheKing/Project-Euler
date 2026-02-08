#include <array>
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

std::array<i64, 10> digit_factorials() {
    std::array<i64, 10> f{};
    f[0] = 1;
    for (int d = 1; d <= 9; ++d) {
        f[static_cast<std::size_t>(d)] = f[static_cast<std::size_t>(d - 1)] * d;
    }
    return f;
}

bool is_curious(const i64 value, const std::array<i64, 10>& fact) {
    i64 x = value;
    i64 sum = 0;
    while (x > 0) {
        sum += fact[static_cast<std::size_t>(x % 10)];
        x /= 10;
    }
    return sum == value;
}

i64 solve() {
    const auto fact = digit_factorials();
    const i64 nine_factorial = fact[9];

    int digits = 1;
    i64 pow10 = 1;
    while (pow10 <= static_cast<i64>(digits) * nine_factorial) {
        ++digits;
        pow10 *= 10;
    }
    const i64 upper = static_cast<i64>(digits - 1) * nine_factorial;

    i64 total = 0;
    for (i64 n = 10; n <= upper; ++n) {
        if (is_curious(n, fact)) {
            total += n;
        }
    }
    return total;
}

bool run_checkpoints() {
    const auto fact = digit_factorials();
    if (!is_curious(145, fact)) {
        std::cerr << "Checkpoint failed for 145" << '\n';
        return false;
    }
    if (!is_curious(40585, fact)) {
        std::cerr << "Checkpoint failed for 40585" << '\n';
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
