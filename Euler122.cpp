#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int limit = 200;
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
    return options.limit >= 1;
}

void dfs_chains(const int max_depth,
                const int limit,
                std::vector<int>& chain,
                std::vector<int>& best) {
    const int depth = static_cast<int>(chain.size()) - 1;
    if (depth >= max_depth) {
        return;
    }

    const int last = chain.back();
    int seen_prev = -1;

    for (int i = depth; i >= 0; --i) {
        const int next = last + chain[static_cast<std::size_t>(i)];
        if (next > limit || next <= last) {
            continue;
        }
        if (next == seen_prev) {
            continue;
        }
        seen_prev = next;

        if (depth + 1 > best[static_cast<std::size_t>(next)]) {
            continue;
        }

        if (depth + 1 < best[static_cast<std::size_t>(next)]) {
            best[static_cast<std::size_t>(next)] = depth + 1;
        }
        chain.push_back(next);
        dfs_chains(max_depth, limit, chain, best);
        chain.pop_back();
    }
}

std::vector<int> minimal_multiplications(const int limit) {
    std::vector<int> best(static_cast<std::size_t>(limit + 1), std::numeric_limits<int>::max());
    best[1] = 0;

    std::vector<int> chain = {1};

    for (int max_depth = 1; max_depth <= 15; ++max_depth) {
        dfs_chains(max_depth, limit, chain, best);

        bool done = true;
        for (int k = 1; k <= limit; ++k) {
            if (best[static_cast<std::size_t>(k)] == std::numeric_limits<int>::max()) {
                done = false;
                break;
            }
        }
        if (done) {
            break;
        }
    }

    return best;
}

int solve(const int limit) {
    const std::vector<int> best = minimal_multiplications(limit);
    int sum = 0;
    for (int k = 1; k <= limit; ++k) {
        sum += best[static_cast<std::size_t>(k)];
    }
    return sum;
}

bool run_checkpoints() {
    const std::vector<int> best = minimal_multiplications(15);
    if (best[15] != 5) {
        std::cerr << "Checkpoint failed for m(15)" << '\n';
        return false;
    }
    if (solve(10) != 26) {
        std::cerr << "Checkpoint failed for sum up to 10" << '\n';
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
