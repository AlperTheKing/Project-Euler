#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string file = "resources/documents/0079_keylog.txt";
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

std::string derive_passcode(const std::vector<std::string>& logs) {
    std::vector<std::set<int>> edges(10);
    std::vector<int> indeg(10, 0);
    std::vector<bool> present(10, false);

    for (const std::string& s : logs) {
        if (s.size() != 3U) {
            continue;
        }

        const int a = s[0] - '0';
        const int b = s[1] - '0';
        const int c = s[2] - '0';
        present[a] = present[b] = present[c] = true;

        if (edges[a].insert(b).second) {
            ++indeg[b];
        }
        if (edges[b].insert(c).second) {
            ++indeg[c];
        }
    }

    std::priority_queue<int, std::vector<int>, std::greater<int>> q;
    for (int d = 0; d <= 9; ++d) {
        if (present[d] && indeg[d] == 0) {
            q.push(d);
        }
    }

    std::string out;
    while (!q.empty()) {
        const int u = q.top();
        q.pop();
        out.push_back(static_cast<char>('0' + u));

        for (const int v : edges[static_cast<std::size_t>(u)]) {
            --indeg[v];
            if (indeg[v] == 0) {
                q.push(v);
            }
        }
    }

    return out;
}

std::string solve(const std::string& file_path) {
    std::ifstream input(file_path);
    if (!input) {
        throw std::runtime_error("Could not open keylog file: " + file_path);
    }

    std::vector<std::string> logs;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) {
            logs.push_back(line);
        }
    }

    return derive_passcode(logs);
}

bool run_checkpoints() {
    const std::vector<std::string> sample = {"123", "135", "245"};
    if (derive_passcode(sample) != "12345") {
        std::cerr << "Checkpoint failed for toy keylog" << '\n';
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
