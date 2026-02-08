#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;

struct Options {
    int n = 50;
    bool run_checkpoints = true;
};

struct Point {
    int x = 0;
    int y = 0;
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
    return options.n >= 1;
}

int solve(const int n) {
    std::vector<Point> points;
    points.reserve(static_cast<std::size_t>(n * n + 2 * n));

    for (int x = 0; x <= n; ++x) {
        for (int y = 0; y <= n; ++y) {
            if (x == 0 && y == 0) {
                continue;
            }
            points.push_back({x, y});
        }
    }

    int count = 0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        for (std::size_t j = i + 1; j < points.size(); ++j) {
            const Point p = points[i];
            const Point q = points[j];

            const i64 dot_o = static_cast<i64>(p.x) * q.x + static_cast<i64>(p.y) * q.y;
            const i64 dot_p = static_cast<i64>(p.x) * (p.x - q.x) + static_cast<i64>(p.y) * (p.y - q.y);
            const i64 dot_q = static_cast<i64>(q.x) * (q.x - p.x) + static_cast<i64>(q.y) * (q.y - p.y);

            if (dot_o == 0 || dot_p == 0 || dot_q == 0) {
                ++count;
            }
        }
    }

    return count;
}

bool run_checkpoints() {
    if (solve(2) != 14) {
        std::cerr << "Checkpoint failed for n=2" << '\n';
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
