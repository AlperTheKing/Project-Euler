#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0081_matrix.txt";
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

std::vector<std::vector<int>> parse_matrix_csv(const std::string& text) {
    std::vector<std::vector<int>> matrix;
    std::istringstream input(text);
    std::string line;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<int> row;
        std::string token;
        std::istringstream row_stream(line);
        while (std::getline(row_stream, token, ',')) {
            row.push_back(std::stoi(token));
        }
        if (!row.empty()) {
            matrix.push_back(row);
        }
    }

    return matrix;
}

long long min_path_sum_two_ways(const std::vector<std::vector<int>>& matrix) {
    const int n = static_cast<int>(matrix.size());
    const int m = static_cast<int>(matrix[0].size());

    std::vector<std::vector<long long>> dp(
        static_cast<std::size_t>(n), std::vector<long long>(static_cast<std::size_t>(m), 0LL));

    dp[0][0] = matrix[0][0];
    for (int i = 1; i < n; ++i) {
        dp[static_cast<std::size_t>(i)][0] = dp[static_cast<std::size_t>(i - 1)][0] + matrix[static_cast<std::size_t>(i)][0];
    }
    for (int j = 1; j < m; ++j) {
        dp[0][static_cast<std::size_t>(j)] = dp[0][static_cast<std::size_t>(j - 1)] + matrix[0][static_cast<std::size_t>(j)];
    }

    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            dp[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
                matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +
                std::min(dp[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)],
                         dp[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)]);
        }
    }

    return dp[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(m - 1)];
}

long long solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open matrix file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return min_path_sum_two_ways(parse_matrix_csv(buffer.str()));
}

bool run_checkpoints() {
    const std::string sample =
        "131,673,234,103,18\n"
        "201,96,342,965,150\n"
        "630,803,746,422,111\n"
        "537,699,497,121,956\n"
        "805,732,524,37,331\n";

    if (min_path_sum_two_ways(parse_matrix_csv(sample)) != 2427LL) {
        std::cerr << "Checkpoint failed for sample matrix" << '\n';
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
