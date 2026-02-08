#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0105_sets.txt";
    bool run_checkpoints = true;
};

bool parse_string_after_prefix(const std::string& arg,
                               const std::string& prefix,
                               std::string& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    value = tail;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_string_after_prefix(arg, "--file=", options.file)) {
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

bool is_special_sum_set(std::vector<int> values) {
    std::sort(values.begin(), values.end());
    if (!satisfies_size_order_rule(values)) {
        return false;
    }
    if (!has_unique_disjoint_subset_sums(values)) {
        return false;
    }
    return true;
}

std::vector<int> parse_set_line(const std::string& line) {
    std::vector<int> values;
    std::istringstream input(line);
    std::string token;
    while (std::getline(input, token, ',')) {
        values.push_back(std::stoi(token));
    }
    return values;
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open sets file: " + file_path);
    }

    int total = 0;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const std::vector<int> values = parse_set_line(line);
        if (is_special_sum_set(values)) {
            for (int v : values) {
                total += v;
            }
        }
    }

    return total;
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

    try {
        std::cout << solve(options.file) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
