#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0067_triangle.txt";
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

std::vector<std::vector<int>> parse_triangle(const std::string& text) {
    std::vector<std::vector<int>> triangle;
    std::istringstream input(text);
    std::string line;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        std::istringstream row_stream(line);
        std::vector<int> row;
        int value = 0;
        while (row_stream >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            triangle.push_back(row);
        }
    }

    return triangle;
}

int max_path_sum(std::vector<std::vector<int>> triangle) {
    if (triangle.empty()) {
        return 0;
    }

    for (int row = static_cast<int>(triangle.size()) - 2; row >= 0; --row) {
        for (int col = 0; col <= row; ++col) {
            triangle[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] +=
                std::max(triangle[static_cast<std::size_t>(row + 1)][static_cast<std::size_t>(col)],
                         triangle[static_cast<std::size_t>(row + 1)][static_cast<std::size_t>(col + 1)]);
        }
    }

    return triangle[0][0];
}

int solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open triangle file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return max_path_sum(parse_triangle(buffer.str()));
}

bool run_checkpoints() {
    const std::string sample = "3\n7 4\n2 4 6\n8 5 9 3\n";
    if (max_path_sum(parse_triangle(sample)) != 23) {
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

    try {
        std::cout << solve(options.file) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 3;
    }

    return 0;
}
