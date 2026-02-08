#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0107_network.txt";
    bool run_checkpoints = true;
};

struct Edge {
    int u = 0;
    int v = 0;
    int w = 0;
};

struct DisjointSet {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit DisjointSet(int n) : parent(static_cast<std::size_t>(n)), rank(static_cast<std::size_t>(n), 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        if (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] = find(parent[static_cast<std::size_t>(x)]);
        }
        return parent[static_cast<std::size_t>(x)];
    }

    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) {
            return false;
        }
        if (rank[static_cast<std::size_t>(a)] < rank[static_cast<std::size_t>(b)]) {
            std::swap(a, b);
        }
        parent[static_cast<std::size_t>(b)] = a;
        if (rank[static_cast<std::size_t>(a)] == rank[static_cast<std::size_t>(b)]) {
            ++rank[static_cast<std::size_t>(a)];
        }
        return true;
    }
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

long long solve_from_text(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    std::vector<std::vector<std::string>> rows;

    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> row;
        std::istringstream row_input(line);
        std::string token;
        while (std::getline(row_input, token, ',')) {
            row.push_back(token);
        }
        rows.push_back(row);
    }

    const int n = static_cast<int>(rows.size());
    std::vector<Edge> edges;
    long long total_weight = 0;

    for (int i = 0; i < n; ++i) {
        if (static_cast<int>(rows[static_cast<std::size_t>(i)].size()) != n) {
            throw std::runtime_error("Malformed network matrix");
        }
        for (int j = i + 1; j < n; ++j) {
            const std::string& token = rows[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            if (token == "-") {
                continue;
            }
            const int w = std::stoi(token);
            total_weight += w;
            edges.push_back({i, j, w});
        }
    }

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) { return a.w < b.w; });

    DisjointSet dsu(n);
    long long mst_weight = 0;
    int used = 0;

    for (const Edge& e : edges) {
        if (dsu.unite(e.u, e.v)) {
            mst_weight += e.w;
            ++used;
            if (used == n - 1) {
                break;
            }
        }
    }

    return total_weight - mst_weight;
}

long long solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open network file: " + file_path);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return solve_from_text(buffer.str());
}

bool run_checkpoints() {
    const std::string sample =
        "-,16,12,21,-,-,-\n"
        "16,-,-,17,20,-,-\n"
        "12,-,-,28,-,31,-\n"
        "21,17,28,-,18,19,23\n"
        "-,20,-,18,-,-,11\n"
        "-,-,31,19,-,-,27\n"
        "-,-,-,23,11,27,-\n";

    if (solve_from_text(sample) != 150) {
        std::cerr << "Checkpoint failed for statement sample network" << '\n';
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
