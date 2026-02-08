#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    int width = 32;
    int height = 10;
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
        if (parse_int_after_prefix(arg, "--width=", options.width) ||
            parse_int_after_prefix(arg, "--height=", options.height)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.width >= 2 && options.height >= 1;
}

std::vector<std::uint32_t> generate_rows(const int width) {
    std::vector<std::uint32_t> rows;

    const auto dfs = [&](auto&& self, int pos, std::uint32_t mask) -> void {
        if (pos == width) {
            rows.push_back(mask);
            return;
        }

        for (int brick : {2, 3}) {
            const int next = pos + brick;
            if (next > width) {
                continue;
            }
            std::uint32_t next_mask = mask;
            if (next < width) {
                next_mask |= (1U << next);
            }
            self(self, next, next_mask);
        }
    };

    dfs(dfs, 0, 0U);
    return rows;
}

u64 solve(const int width, const int height) {
    const std::vector<std::uint32_t> rows = generate_rows(width);
    const int n = static_cast<int>(rows.size());

    std::vector<std::vector<int>> compatible(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if ((rows[static_cast<std::size_t>(i)] & rows[static_cast<std::size_t>(j)]) == 0U) {
                compatible[static_cast<std::size_t>(i)].push_back(j);
            }
        }
    }

    std::vector<u64> dp(static_cast<std::size_t>(n), 1ULL);

    for (int h = 1; h < height; ++h) {
        std::vector<u64> next(static_cast<std::size_t>(n), 0ULL);
        for (int i = 0; i < n; ++i) {
            u64 ways = 0;
            for (int j : compatible[static_cast<std::size_t>(i)]) {
                ways += dp[static_cast<std::size_t>(j)];
            }
            next[static_cast<std::size_t>(i)] = ways;
        }
        dp.swap(next);
    }

    u64 total = 0;
    for (u64 v : dp) {
        total += v;
    }
    return total;
}

bool run_checkpoints() {
    if (solve(9, 3) != 8ULL) {
        std::cerr << "Checkpoint failed for width=9,height=3" << '\n';
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

    std::cout << solve(options.width, options.height) << '\n';
    return 0;
}
