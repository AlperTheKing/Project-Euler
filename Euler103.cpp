#include <algorithm>
#include <iostream>
#include <numeric>
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

bool satisfies_size_order_rule(const std::vector<int>& values) {
    const int n = static_cast<int>(values.size());
    std::vector<int> prefix(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        prefix[static_cast<std::size_t>(i + 1)] = prefix[static_cast<std::size_t>(i)] + values[static_cast<std::size_t>(i)];
    }

    for (int m = 1; 2 * m + 1 <= n; ++m) {
        const int left = prefix[static_cast<std::size_t>(m + 1)] - prefix[0];
        const int right = prefix[static_cast<std::size_t>(n)] - prefix[static_cast<std::size_t>(n - m)];
        if (left <= right) {
            return false;
        }
    }

    return true;
}

bool has_unique_disjoint_subset_sums(const std::vector<int>& values) {
    const int n = static_cast<int>(values.size());
    const int total_masks = 1 << n;
    std::vector<int> subset_sum(static_cast<std::size_t>(total_masks), 0);

    for (int mask = 1; mask < total_masks; ++mask) {
        const int bit = __builtin_ctz(static_cast<unsigned>(mask));
        const int prev = mask & (mask - 1);
        subset_sum[static_cast<std::size_t>(mask)] =
            subset_sum[static_cast<std::size_t>(prev)] + values[static_cast<std::size_t>(bit)];
    }

    for (int a = 1; a < total_masks; ++a) {
        for (int b = a + 1; b < total_masks; ++b) {
            if ((a & b) != 0) {
                continue;
            }
            if (subset_sum[static_cast<std::size_t>(a)] == subset_sum[static_cast<std::size_t>(b)]) {
                return false;
            }
        }
    }

    return true;
}

bool is_special_sum_set(const std::vector<int>& values) {
    if (!std::is_sorted(values.begin(), values.end())) {
        return false;
    }
    if (!satisfies_size_order_rule(values)) {
        return false;
    }
    if (!has_unique_disjoint_subset_sums(values)) {
        return false;
    }
    return true;
}

bool partial_size_order_prune(const std::vector<int>& values) {
    const int n = static_cast<int>(values.size());
    std::vector<int> prefix(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        prefix[static_cast<std::size_t>(i + 1)] = prefix[static_cast<std::size_t>(i)] + values[static_cast<std::size_t>(i)];
    }

    for (int m = 1; 2 * m + 1 <= n; ++m) {
        const int left = prefix[static_cast<std::size_t>(m + 1)] - prefix[0];
        const int right = prefix[static_cast<std::size_t>(n)] - prefix[static_cast<std::size_t>(n - m)];
        if (left <= right) {
            return false;
        }
    }

    return true;
}

std::vector<int> near_optimum_seed_for_seven() {
    const std::vector<int> previous_optimum = {11, 18, 19, 20, 22, 25};
    const int middle = previous_optimum[previous_optimum.size() / 2];

    std::vector<int> seed;
    seed.push_back(middle);
    for (int v : previous_optimum) {
        seed.push_back(middle + v);
    }
    return seed;
}

void search_best_set(const int target_size,
                     int index,
                     int next_value,
                     int current_sum,
                     std::vector<int>& current,
                     std::vector<int>& best,
                     int& best_sum) {
    const int remaining = target_size - index;
    const int min_possible_tail = remaining * next_value + (remaining * (remaining - 1)) / 2;
    if (current_sum + min_possible_tail >= best_sum) {
        return;
    }

    if (index == target_size) {
        if (is_special_sum_set(current)) {
            best_sum = current_sum;
            best = current;
        }
        return;
    }

    for (int v = next_value;; ++v) {
        const int rest = target_size - index - 1;
        const int min_rest = rest * (v + 1) + (rest * (rest - 1)) / 2;
        if (current_sum + v + min_rest >= best_sum) {
            break;
        }

        current.push_back(v);
        if (partial_size_order_prune(current)) {
            search_best_set(target_size, index + 1, v + 1, current_sum + v, current, best, best_sum);
        }
        current.pop_back();
    }
}

std::string solve() {
    std::vector<int> best = near_optimum_seed_for_seven();
    int best_sum = std::accumulate(best.begin(), best.end(), 0);

    std::vector<int> current;
    search_best_set(7, 0, 1, 0, current, best, best_sum);

    std::string concat;
    for (int v : best) {
        concat += std::to_string(v);
    }
    return concat;
}

bool run_checkpoints() {
    if (!is_special_sum_set({3, 5, 6, 7})) {
        std::cerr << "Checkpoint failed: known special set rejected" << '\n';
        return false;
    }
    if (is_special_sum_set({2, 3, 4, 5})) {
        std::cerr << "Checkpoint failed: equal-sum disjoint subsets not detected" << '\n';
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
