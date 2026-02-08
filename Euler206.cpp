#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

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

bool matches_pattern(u64 sq) {
    if (sq % 10ULL != 0ULL) {
        return false;
    }
    sq /= 100ULL;

    for (int d = 9; d >= 1; --d) {
        if (sq % 10ULL != static_cast<u64>(d)) {
            return false;
        }
        sq /= 100ULL;
    }

    return true;
}

i64 solve() {
    const u64 low_sq = 1020304050607080900ULL;
    const u64 high_sq = 1929394959697989990ULL;

    i64 low = static_cast<i64>(std::ceil(std::sqrt(static_cast<long double>(low_sq))));
    const i64 high = static_cast<i64>(std::floor(std::sqrt(static_cast<long double>(high_sq))));

    while (low % 10 != 0) {
        ++low;
    }

    for (i64 base = low; base <= high; base += 10) {
        const int tail = static_cast<int>(base % 100);
        if (tail != 30 && tail != 70) {
            continue;
        }

        const u64 sq = static_cast<u64>(base) * static_cast<u64>(base);
        if (matches_pattern(sq)) {
            return base;
        }
    }

    return -1;
}

bool run_checkpoints() {
    if (!matches_pattern(1020304050607080900ULL)) {
        std::cerr << "Checkpoint failed for pattern matcher positive case" << '\n';
        return false;
    }
    if (matches_pattern(1020304050607080000ULL)) {
        std::cerr << "Checkpoint failed for pattern matcher negative case" << '\n';
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
