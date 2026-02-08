#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int limit = 1000000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.limit >= 2;
}

std::vector<int> proper_divisor_sums(const int limit) {
    std::vector<int> s(static_cast<std::size_t>(limit + 1), 0);
    for (int d = 1; d <= limit / 2; ++d) {
        for (int m = 2 * d; m <= limit; m += d) {
            s[static_cast<std::size_t>(m)] += d;
        }
    }
    return s;
}

int solve(const int limit) {
    const std::vector<int> sum = proper_divisor_sums(limit);

    std::vector<char> done(static_cast<std::size_t>(limit + 1), 0);
    int best_len = 0;
    int best_min_member = 0;

    for (int start = 2; start <= limit; ++start) {
        if (done[static_cast<std::size_t>(start)]) {
            continue;
        }

        std::vector<int> path;
        std::unordered_map<int, int> index;
        int n = start;

        while (n >= 1 && n <= limit && index.find(n) == index.end() && !done[static_cast<std::size_t>(n)]) {
            index[n] = static_cast<int>(path.size());
            path.push_back(n);
            n = sum[static_cast<std::size_t>(n)];
        }

        if (n >= 1 && n <= limit) {
            auto it = index.find(n);
            if (it != index.end()) {
                const int loop_start = it->second;
                const int loop_len = static_cast<int>(path.size()) - loop_start;
                int min_member = path[static_cast<std::size_t>(loop_start)];
                for (int i = loop_start; i < static_cast<int>(path.size()); ++i) {
                    min_member = std::min(min_member, path[static_cast<std::size_t>(i)]);
                }

                if (loop_len > best_len || (loop_len == best_len && min_member < best_min_member)) {
                    best_len = loop_len;
                    best_min_member = min_member;
                }
            }
        }

        for (const int v : path) {
            done[static_cast<std::size_t>(v)] = 1;
        }
    }

    return best_min_member;
}

bool run_checkpoints() {
    if (solve(300) != 220) {
        std::cerr << "Checkpoint failed for limit=300" << '\n';
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

    std::cout << solve(options.limit) << '\n';
    return 0;
}
