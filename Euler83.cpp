#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0083_matrix.txt";
    bool run_checkpoints = true;
};

struct Node {
    int r = 0;
    int c = 0;
    long long dist = 0;
};

struct NodeCmp {
    bool operator()(const Node& a, const Node& b) const { return a.dist > b.dist; }
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

long long min_path_sum_four_ways(const std::vector<std::vector<int>>& matrix) {
    const int n = static_cast<int>(matrix.size());
    const int m = static_cast<int>(matrix[0].size());

    std::vector<std::vector<long long>> dist(
        static_cast<std::size_t>(n),
        std::vector<long long>(static_cast<std::size_t>(m), std::numeric_limits<long long>::max()));

    std::priority_queue<Node, std::vector<Node>, NodeCmp> pq;
    dist[0][0] = matrix[0][0];
    pq.push({0, 0, dist[0][0]});

    const int dr[4] = {1, -1, 0, 0};
    const int dc[4] = {0, 0, 1, -1};

    while (!pq.empty()) {
        const Node cur = pq.top();
        pq.pop();

        if (cur.dist != dist[static_cast<std::size_t>(cur.r)][static_cast<std::size_t>(cur.c)]) {
            continue;
        }
        if (cur.r == n - 1 && cur.c == m - 1) {
            return cur.dist;
        }

        for (int dir = 0; dir < 4; ++dir) {
            const int nr = cur.r + dr[dir];
            const int nc = cur.c + dc[dir];
            if (nr < 0 || nr >= n || nc < 0 || nc >= m) {
                continue;
            }

            const long long nd = cur.dist + matrix[static_cast<std::size_t>(nr)][static_cast<std::size_t>(nc)];
            if (nd < dist[static_cast<std::size_t>(nr)][static_cast<std::size_t>(nc)]) {
                dist[static_cast<std::size_t>(nr)][static_cast<std::size_t>(nc)] = nd;
                pq.push({nr, nc, nd});
            }
        }
    }

    return dist[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(m - 1)];
}

long long solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open matrix file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return min_path_sum_four_ways(parse_matrix_csv(buffer.str()));
}

bool run_checkpoints() {
    const std::string sample =
        "131,673,234,103,18\n"
        "201,96,342,965,150\n"
        "630,803,746,422,111\n"
        "537,699,497,121,956\n"
        "805,732,524,37,331\n";

    if (min_path_sum_four_ways(parse_matrix_csv(sample)) != 2297LL) {
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
