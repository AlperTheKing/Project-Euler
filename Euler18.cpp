#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    bool run_checkpoints = true;
};

const std::vector<std::vector<int>> kTriangle = {
    {75},
    {95, 64},
    {17, 47, 82},
    {18, 35, 87, 10},
    {20, 4, 82, 47, 65},
    {19, 1, 23, 75, 3, 34},
    {88, 2, 77, 73, 7, 63, 67},
    {99, 65, 4, 28, 6, 16, 70, 92},
    {41, 41, 26, 56, 83, 40, 80, 70, 33},
    {41, 48, 72, 33, 47, 32, 37, 16, 94, 29},
    {53, 71, 44, 65, 25, 43, 91, 52, 97, 51, 14},
    {70, 11, 33, 28, 77, 73, 17, 78, 39, 68, 17, 57},
    {91, 71, 52, 38, 17, 14, 91, 43, 58, 50, 27, 29, 48},
    {63, 66, 4, 68, 89, 53, 67, 30, 73, 16, 69, 87, 40, 31},
    {4, 62, 98, 27, 23, 9, 70, 98, 73, 93, 38, 53, 60, 4, 23},
};

bool parse_arguments(const int argc, char** argv, Options& options) {
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

int max_path_sum(const std::vector<std::vector<int>>& triangle) {
    if (triangle.empty()) {
        return 0;
    }

    std::vector<int> dp = triangle.back();
    for (int row = static_cast<int>(triangle.size()) - 2; row >= 0; --row) {
        for (int col = 0; col <= row; ++col) {
            dp[col] = triangle[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] +
                      std::max(dp[col], dp[col + 1]);
        }
    }

    return dp[0];
}

int solve() {
    return max_path_sum(kTriangle);
}

bool run_checkpoints() {
    const std::vector<std::vector<int>> sample = {
        {3},
        {7, 4},
        {2, 4, 6},
        {8, 5, 9, 3},
    };

    if (max_path_sum(sample) != 23) {
        std::cerr << "Checkpoint failed for sample triangle" << '\n';
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
