#include <algorithm>
#include <array>
#include <iostream>
#include <set>
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

unsigned long long solve() {
    std::array<int, 10> nums{};
    for (int i = 0; i < 10; ++i) {
        nums[static_cast<std::size_t>(i)] = i + 1;
    }

    std::string best;

    // Choose inner ring (5 numbers in cyclic order), derive outer ring from line sums.
    std::array<int, 10> mask{};
    std::fill(mask.begin(), mask.begin() + 5, 1);

    do {
        std::vector<int> inner_set;
        std::vector<int> outer_set;
        inner_set.reserve(5);
        outer_set.reserve(5);

        for (int i = 0; i < 10; ++i) {
            if (mask[static_cast<std::size_t>(i)] == 1) {
                inner_set.push_back(nums[static_cast<std::size_t>(i)]);
            } else {
                outer_set.push_back(nums[static_cast<std::size_t>(i)]);
            }
        }

        std::sort(inner_set.begin(), inner_set.end());
        do {
            const int i1 = inner_set[0];
            const int i2 = inner_set[1];
            const int i3 = inner_set[2];
            const int i4 = inner_set[3];
            const int i5 = inner_set[4];

            for (const int o1 : outer_set) {
                const int target_sum = o1 + i1 + i2;

                const int o2 = target_sum - i2 - i3;
                const int o3 = target_sum - i3 - i4;
                const int o4 = target_sum - i4 - i5;
                const int o5 = target_sum - i5 - i1;

                std::array<int, 5> outer = {o1, o2, o3, o4, o5};
                std::set<int> outer_needed(outer_set.begin(), outer_set.end());

                bool ok = true;
                for (const int o : outer) {
                    if (outer_needed.erase(o) == 0U) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) {
                    continue;
                }

                const int min_outer = *std::min_element(outer.begin(), outer.end());
                if (outer[0] != min_outer) {
                    continue;
                }

                std::string repr;
                repr.reserve(20U);
                auto append_line = [&](int a, int b, int c) {
                    repr += std::to_string(a);
                    repr += std::to_string(b);
                    repr += std::to_string(c);
                };

                append_line(outer[0], i1, i2);
                append_line(outer[1], i2, i3);
                append_line(outer[2], i3, i4);
                append_line(outer[3], i4, i5);
                append_line(outer[4], i5, i1);

                if (repr.size() == 16U && repr > best) {
                    best = repr;
                }
            }
        } while (std::next_permutation(inner_set.begin(), inner_set.end()));

    } while (std::prev_permutation(mask.begin(), mask.end()));

    return std::stoull(best);
}

bool run_checkpoints() {
    const unsigned long long value = solve();
    if (value == 0ULL) {
        std::cerr << "Checkpoint failed for non-zero result" << '\n';
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
