#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

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

i64 polygonal(const int type, const int n) {
    switch (type) {
        case 3:
            return static_cast<i64>(n) * (n + 1) / 2;
        case 4:
            return static_cast<i64>(n) * n;
        case 5:
            return static_cast<i64>(n) * (3 * n - 1) / 2;
        case 6:
            return static_cast<i64>(n) * (2 * n - 1);
        case 7:
            return static_cast<i64>(n) * (5 * n - 3) / 2;
        case 8:
            return static_cast<i64>(n) * (3 * n - 2);
        default:
            return -1;
    }
}

std::vector<int> generate_four_digit_values(const int type) {
    std::vector<int> values;
    for (int n = 1;; ++n) {
        const i64 x = polygonal(type, n);
        if (x > 9999) {
            break;
        }
        if (x >= 1000) {
            const int value = static_cast<int>(x);
            const int prefix = value / 100;
            const int suffix = value % 100;
            if (prefix >= 10 && suffix >= 10) {
                values.push_back(value);
            }
        }
    }
    return values;
}

int solve() {
    const std::array<int, 6> types = {3, 4, 5, 6, 7, 8};
    std::array<std::vector<int>, 6> nums;
    for (int i = 0; i < 6; ++i) {
        nums[static_cast<std::size_t>(i)] = generate_four_digit_values(types[static_cast<std::size_t>(i)]);
    }

    std::vector<int> chosen_values;
    std::vector<int> chosen_types;
    int answer = -1;

    auto dfs = [&](auto&& self, int suffix_needed, int used_mask, int first_prefix, int current_sum) -> void {
        if (answer != -1) {
            return;
        }
        if (static_cast<int>(chosen_values.size()) == 6) {
            if (suffix_needed == first_prefix) {
                answer = current_sum;
            }
            return;
        }

        for (int t = 0; t < 6; ++t) {
            if ((used_mask >> t) & 1) {
                continue;
            }

            for (const int value : nums[static_cast<std::size_t>(t)]) {
                const int prefix = value / 100;
                const int suffix = value % 100;

                if (!chosen_values.empty() && prefix != suffix_needed) {
                    continue;
                }

                chosen_values.push_back(value);
                chosen_types.push_back(t);

                const int next_first_prefix = chosen_values.size() == 1U ? prefix : first_prefix;
                self(self, suffix, used_mask | (1 << t), next_first_prefix, current_sum + value);

                chosen_types.pop_back();
                chosen_values.pop_back();
            }
        }
    };

    dfs(dfs, -1, 0, -1, 0);
    return answer;
}

bool run_checkpoints() {
    if (polygonal(3, 45) != 1035 || polygonal(8, 19) != 1045) {
        std::cerr << "Checkpoint failed for polygonal formulas" << '\n';
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
