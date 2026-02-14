#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <functional>

namespace {

struct Options {
    int players = 100;
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
        if (parse_int_after_prefix(arg, "--players=", options.players)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.players >= 2 && (options.players % 2 == 0);
}

long double solve(const int players) {
    const std::array<int, 5> delta = {-2, -1, 0, 1, 2};
    const std::array<long double, 5> prob = {
        1.0L / 36.0L,
        2.0L / 9.0L,
        1.0L / 2.0L,
        2.0L / 9.0L,
        1.0L / 36.0L,
    };

    const int n = players;
    const int m = n - 1;

    std::vector<std::vector<long double>> a(static_cast<std::size_t>(m),
                                            std::vector<long double>(static_cast<std::size_t>(m + 1), 0.0L));

    for (int k = 1; k < n; ++k) {
        const int row = k - 1;
        a[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)] = 1.0L;
        a[static_cast<std::size_t>(row)][static_cast<std::size_t>(m)] = 1.0L;

        for (int i = 0; i < 5; ++i) {
            int nk = (k + delta[static_cast<std::size_t>(i)]) % n;
            if (nk < 0) {
                nk += n;
            }
            if (nk == 0) {
                continue;
            }
            a[static_cast<std::size_t>(row)][static_cast<std::size_t>(nk - 1)] -= prob[static_cast<std::size_t>(i)];
        }
    }

    for (int col = 0; col < m; ++col) {
        int pivot = col;
        for (int row = col + 1; row < m; ++row) {
            if (std::fabsl(a[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]) >
                std::fabsl(a[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(col)])) {
                pivot = row;
            }
        }
        std::swap(a[static_cast<std::size_t>(col)], a[static_cast<std::size_t>(pivot)]);

        const long double div = a[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
        for (int j = col; j <= m; ++j) {
            a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)] /= div;
        }

        for (int row = 0; row < m; ++row) {
            if (row == col) {
                continue;
            }
            const long double factor = a[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
            if (std::fabsl(factor) < 1e-21L) {
                continue;
            }
            for (int j = col; j <= m; ++j) {
                a[static_cast<std::size_t>(row)][static_cast<std::size_t>(j)] -=
                    factor * a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)];
            }
        }
    }

    const int start = n / 2;
    return a[static_cast<std::size_t>(start - 1)][static_cast<std::size_t>(m)];
}

bool run_checkpoints() {
    if (std::fabsl(solve(2) - 2.25L) > 1e-12L) {
        std::cerr << "Checkpoint failed for 2 players" << '\n';
        return false;
    }
    if (std::fabsl(solve(4) - 7.2L) > 1e-12L) {
        std::cerr << "Checkpoint failed for 4 players" << '\n';
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

    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(solve(options.players)) << '\n';
    return 0;
}
