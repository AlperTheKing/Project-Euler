#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = long long;

struct Options {
    int width = 47;
    int height = 43;
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
    return options.width >= 1 && options.height >= 1;
}

i64 solve(const int width, const int height) {
    int n = std::min(width, height);
    int m = std::max(width, height);

    i64 answer = static_cast<i64>(n) * (n + 1) * (n + 2) / 6 *
                 static_cast<i64>(m) * (m + 1) * (m + 2) / 6;

    const int shift = m + 3;
    const int dim = 2 * shift + 2;
    std::vector<std::vector<int>> grid(static_cast<std::size_t>(dim), std::vector<int>(static_cast<std::size_t>(dim), 0));

    for (int i = 0; i <= 2 * n; ++i) {
        for (int j = 0; j <= 2 * m; ++j) {
            if (((i + j) & 1) != 0) {
                continue;
            }
            const int x = (i + j) / 2;
            const int y = (i + 2 * shift - j) / 2;
            grid[static_cast<std::size_t>(x + 1)][static_cast<std::size_t>(y + 1)] = 1;
        }
    }

    for (int i = 1; i < dim; ++i) {
        for (int j = 1; j < dim; ++j) {
            grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
                grid[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j)] +
                grid[static_cast<std::size_t>(i)][static_cast<std::size_t>(j - 1)] -
                grid[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(j - 1)];
        }
    }

    for (int lx = 1; lx < dim; ++lx) {
        for (int rx = lx + 1; rx < dim; ++rx) {
            for (int ly = 1; ly < dim; ++ly) {
                for (int ry = ly + 1; ry < dim; ++ry) {
                    const int area = (rx - lx + 1) * (ry - ly + 1);
                    const int full = grid[static_cast<std::size_t>(rx)][static_cast<std::size_t>(ry)] -
                                     grid[static_cast<std::size_t>(rx)][static_cast<std::size_t>(ly - 1)] -
                                     grid[static_cast<std::size_t>(lx - 1)][static_cast<std::size_t>(ry)] +
                                     grid[static_cast<std::size_t>(lx - 1)][static_cast<std::size_t>(ly - 1)];
                    if (full != area) {
                        continue;
                    }

                    const int mi = (rx + ry - shift - 1) / 2;
                    const int mj = (rx - ly + shift + 1) / 2;
                    if (mi < 1 || mj < 1 || mi > n || mj > m) {
                        continue;
                    }
                    answer += static_cast<i64>(n - mi + 1) * (m - mj + 1);
                }
            }
        }
    }

    return answer;
}

bool run_checkpoints() {
    if (solve(3, 2) != 72LL) {
        std::cerr << "Checkpoint failed for 3x2" << '\n';
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
