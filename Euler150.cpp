#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int rows = 1000;
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
        if (parse_int_after_prefix(arg, "--rows=", options.rows)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.rows >= 1;
}

std::vector<std::vector<int>> generate_triangle(const int rows) {
    std::vector<std::vector<int>> triangle(static_cast<std::size_t>(rows));

    std::int64_t t = 0;
    for (int r = 0; r < rows; ++r) {
        triangle[static_cast<std::size_t>(r)].reserve(static_cast<std::size_t>(r + 1));
        for (int c = 0; c <= r; ++c) {
            t = (615949LL * t + 797807LL) % (1LL << 20);
            triangle[static_cast<std::size_t>(r)].push_back(static_cast<int>(t - (1LL << 19)));
        }
    }

    return triangle;
}

i64 min_subtriangle_sum(const std::vector<std::vector<int>>& triangle) {
    const int rows = static_cast<int>(triangle.size());

    std::vector<std::vector<i64>> prefix(static_cast<std::size_t>(rows));
    for (int r = 0; r < rows; ++r) {
        prefix[static_cast<std::size_t>(r)].resize(static_cast<std::size_t>(r + 2), 0);
        for (int c = 0; c <= r; ++c) {
            prefix[static_cast<std::size_t>(r)][static_cast<std::size_t>(c + 1)] =
                prefix[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] +
                triangle[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
        }
    }

    i64 best = std::numeric_limits<i64>::max();

    for (int top = 0; top < rows; ++top) {
        std::vector<i64> accum(static_cast<std::size_t>(top + 1), 0);

        for (int bottom = top; bottom < rows; ++bottom) {
            const int height = bottom - top;

            for (int c = 0; c <= top; ++c) {
                const i64 row_segment =
                    prefix[static_cast<std::size_t>(bottom)][static_cast<std::size_t>(c + height + 1)] -
                    prefix[static_cast<std::size_t>(bottom)][static_cast<std::size_t>(c)];
                accum[static_cast<std::size_t>(c)] += row_segment;
                best = std::min(best, accum[static_cast<std::size_t>(c)]);
            }
        }
    }

    return best;
}

i64 solve(const int rows) {
    return min_subtriangle_sum(generate_triangle(rows));
}

bool run_checkpoints() {
    const std::vector<std::vector<int>> sample = {
        {15},
        {-14, -7},
        {20, -13, -5},
        {-3, 8, 23, -26},
        {1, -4, -5, -18, 5},
        {-16, 31, 2, 9, 28, 3},
    };

    if (min_subtriangle_sum(sample) != -42) {
        std::cerr << "Checkpoint failed for statement sample triangle" << '\n';
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

    std::cout << solve(options.rows) << '\n';
    return 0;
}
