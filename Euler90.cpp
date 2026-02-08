#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace {

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

bool has_digit(const int mask, const int d) {
    if (d == 6 || d == 9) {
        return ((mask >> 6) & 1) || ((mask >> 9) & 1);
    }
    return ((mask >> d) & 1) != 0;
}

bool can_display_all_squares(const int a_mask, const int b_mask) {
    static const std::array<std::pair<int, int>, 9> squares = {
        std::pair<int, int>{0, 1}, {0, 4}, {0, 9}, {1, 6}, {2, 5},
        {3, 6},                {4, 9}, {6, 4}, {8, 1},
    };

    for (const auto& [x, y] : squares) {
        const bool ok = (has_digit(a_mask, x) && has_digit(b_mask, y)) ||
                        (has_digit(a_mask, y) && has_digit(b_mask, x));
        if (!ok) {
            return false;
        }
    }
    return true;
}

int solve() {
    std::vector<int> masks;
    for (int mask = 0; mask < (1 << 10); ++mask) {
        if (__builtin_popcount(static_cast<unsigned>(mask)) == 6) {
            masks.push_back(mask);
        }
    }

    int count = 0;
    for (std::size_t i = 0; i < masks.size(); ++i) {
        for (std::size_t j = i; j < masks.size(); ++j) {
            if (can_display_all_squares(masks[i], masks[j])) {
                ++count;
            }
        }
    }

    return count;
}

bool run_checkpoints() {
    const int mask_with_6 = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 6);
    if (!has_digit(mask_with_6, 9)) {
        std::cerr << "Checkpoint failed for 6/9 interchange" << '\n';
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
