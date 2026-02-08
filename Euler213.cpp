#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Options {
    int size = 30;
    int steps = 50;
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
        if (parse_int_after_prefix(arg, "--size=", options.size) ||
            parse_int_after_prefix(arg, "--steps=", options.steps)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.size >= 1 && options.steps >= 0;
}

long double solve(const int size, const int steps) {
    const int cells = size * size;

    std::vector<std::vector<int>> neighbors(static_cast<std::size_t>(cells));
    auto index = [&](const int r, const int c) { return r * size + c; };

    for (int r = 0; r < size; ++r) {
        for (int c = 0; c < size; ++c) {
            const int id = index(r, c);
            if (r > 0) {
                neighbors[static_cast<std::size_t>(id)].push_back(index(r - 1, c));
            }
            if (r + 1 < size) {
                neighbors[static_cast<std::size_t>(id)].push_back(index(r + 1, c));
            }
            if (c > 0) {
                neighbors[static_cast<std::size_t>(id)].push_back(index(r, c - 1));
            }
            if (c + 1 < size) {
                neighbors[static_cast<std::size_t>(id)].push_back(index(r, c + 1));
            }
        }
    }

    std::vector<long double> empty_prob(static_cast<std::size_t>(cells), 1.0L);

    std::vector<long double> dist(static_cast<std::size_t>(cells), 0.0L);
    std::vector<long double> next(static_cast<std::size_t>(cells), 0.0L);

    for (int start = 0; start < cells; ++start) {
        std::fill(dist.begin(), dist.end(), 0.0L);
        dist[static_cast<std::size_t>(start)] = 1.0L;

        for (int step = 0; step < steps; ++step) {
            std::fill(next.begin(), next.end(), 0.0L);

            for (int i = 0; i < cells; ++i) {
                const long double p = dist[static_cast<std::size_t>(i)];
                if (p == 0.0L) {
                    continue;
                }
                const long double share = p / static_cast<long double>(neighbors[static_cast<std::size_t>(i)].size());
                for (int nb : neighbors[static_cast<std::size_t>(i)]) {
                    next[static_cast<std::size_t>(nb)] += share;
                }
            }

            dist.swap(next);
        }

        for (int i = 0; i < cells; ++i) {
            empty_prob[static_cast<std::size_t>(i)] *= (1.0L - dist[static_cast<std::size_t>(i)]);
        }
    }

    long double expected = 0.0L;
    for (long double p : empty_prob) {
        expected += p;
    }
    return expected;
}

bool run_checkpoints() {
    if (std::abs(solve(2, 1) - 1.0L) > 1e-15L) {
        std::cerr << "Checkpoint failed for 2x2 after 1 step" << '\n';
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

    std::cout << std::fixed << std::setprecision(6)
              << static_cast<double>(solve(options.size, options.steps)) << '\n';
    return 0;
}
