#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int n = 5;
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
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 2 && options.n <= 10;
}

u64 solve(const int n) {
    const int edge_count = 1 << n;
    const int node_mask = (1 << (n - 1)) - 1;

    std::vector<std::uint8_t> used_edge(static_cast<std::size_t>(edge_count), 0U);
    std::vector<int> sequence(static_cast<std::size_t>(n), 0);
    used_edge[0] = 1U;  // The starting N-bit block is all zeros.

    u64 total = 0;

    const auto dfs = [&](auto&& self, int node, int used) -> void {
        if (used == edge_count) {
            if (node != 0) {
                return;
            }
            u64 value = 0;
            for (int i = 0; i < edge_count; ++i) {
                value = (value << 1) | static_cast<u64>(sequence[static_cast<std::size_t>(i)]);
            }
            total += value;
            return;
        }

        for (int bit = 0; bit <= 1; ++bit) {
            const int edge = ((node << 1) | bit) & (edge_count - 1);
            if (used_edge[static_cast<std::size_t>(edge)] != 0U) {
                continue;
            }
            used_edge[static_cast<std::size_t>(edge)] = 1U;
            sequence.push_back(bit);
            self(self, edge & node_mask, used + 1);
            sequence.pop_back();
            used_edge[static_cast<std::size_t>(edge)] = 0U;
        }
    };
    dfs(dfs, 0, 1);
    return total;
}

bool run_checkpoints() {
    if (solve(3) != 52ULL) {
        std::cerr << "Checkpoint failed for S(3)=52 sample" << '\n';
        return false;
    }
    if (solve(2) != 3ULL) {
        std::cerr << "Checkpoint failed for S(2)=3" << '\n';
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
    std::cout << solve(options.n) << '\n';
    return 0;
}
