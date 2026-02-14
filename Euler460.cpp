#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Options {
    int d = 10'000;
    int window = 30;
    double radial_eps = 0.5;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) return false;
    try {
        out = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_double_after_prefix(const std::string& arg, const std::string& prefix, double& out) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) return false;
    try {
        out = std::stod(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--d=", options.d)) continue;
        if (parse_int_after_prefix(arg, "--window=", options.window)) continue;
        if (parse_double_after_prefix(arg, "--radial-eps=", options.radial_eps)) continue;

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.d >= 2 && (options.d % 2 == 0) && options.window >= 1 && options.radial_eps > 0.0;
}

double solve_fast(const int d, const int window, const double radial_eps) {
    const int w = d / 2;
    const int h = w;
    const double inf = std::numeric_limits<double>::infinity();

    std::vector<std::vector<double>> g(static_cast<std::size_t>(h + 1),
                                       std::vector<double>(static_cast<std::size_t>(w + 1), inf));
    g[1][0] = 0.0;

    std::vector<double> log_y(static_cast<std::size_t>(h + 1), 0.0);
    for (int y = 1; y <= h; ++y) {
        log_y[static_cast<std::size_t>(y)] = std::log(static_cast<double>(y));
    }

    std::vector<int> x_ref(static_cast<std::size_t>(h + 1), 0);
    for (int y = 1; y <= h; ++y) {
        const int inner = w * w - y * y;
        x_ref[static_cast<std::size_t>(y)] = w - static_cast<int>(std::sqrt(static_cast<double>(inner)));
    }

    std::vector<std::vector<int>> source_y(static_cast<std::size_t>(w));
    const double outer2 = (static_cast<double>(w) + radial_eps) * (static_cast<double>(w) + radial_eps);
    const double inner_r = std::max(0.0, static_cast<double>(w) - radial_eps);
    const double inner2 = inner_r * inner_r;
    for (int x0 = 0; x0 < w; ++x0) {
        const double dx = static_cast<double>(w - x0);
        const double dx2 = dx * dx;
        const double low = inner2 - dx2;
        const double high = outer2 - dx2;
        if (high < 1.0) continue;

        const int y_lo = std::max(1, static_cast<int>(std::ceil(std::sqrt(std::max(0.0, low)))));
        const int y_hi = std::min(h, static_cast<int>(std::floor(std::sqrt(std::max(0.0, high)))));
        if (y_lo > y_hi) continue;

        auto& ys = source_y[static_cast<std::size_t>(x0)];
        ys.reserve(static_cast<std::size_t>(y_hi - y_lo + 1));
        for (int y0 = y_lo; y0 <= y_hi; ++y0) {
            const double r = std::sqrt(dx2 + static_cast<double>(y0) * static_cast<double>(y0));
            if (std::fabs(r - static_cast<double>(w)) <= radial_eps + 1e-12) {
                ys.push_back(y0);
            }
        }
    }

    double best_half = inf;

    for (int x0 = 0; x0 < w; ++x0) {
        const auto& ys = source_y[static_cast<std::size_t>(x0)];
        for (const int y0 : ys) {
            const double base = g[static_cast<std::size_t>(y0)][static_cast<std::size_t>(x0)];
            if (!std::isfinite(base)) continue;

            for (int y = y0; y <= h; ++y) {
                const int dy = y - y0;
                const double inv_v = (dy == 0)
                                         ? 1.0 / static_cast<double>(y0)
                                         : (log_y[static_cast<std::size_t>(y)] -
                                            log_y[static_cast<std::size_t>(y0)]) /
                                               static_cast<double>(dy);

                const int xr = x_ref[static_cast<std::size_t>(y)];
                const int lx = std::max(x0, xr - window);
                const int rx = std::min(w, xr + 1);

                for (int x = lx; x <= rx; ++x) {
                    const int dx = x - x0;
                    const double len = std::sqrt(static_cast<double>(dx * dx + dy * dy)) * inv_v;
                    double& dst = g[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)];
                    const double cand = base + len;
                    if (cand < dst) {
                        dst = cand;
                        if (x == w && cand < best_half) best_half = cand;
                    }
                }
            }
        }
    }

    return 2.0 * best_half;
}

bool close_to(const double a, const double b, const double eps) {
    return std::fabs(a - b) <= eps;
}

bool run_checkpoints(const int window, const double radial_eps) {
    if (!close_to(solve_fast(4, window, radial_eps), 2.960516287, 5e-10)) {
        std::cerr << "Checkpoint failed: F(4)\n";
        return false;
    }
    if (!close_to(solve_fast(10, window, radial_eps), 4.668187834, 5e-10)) {
        std::cerr << "Checkpoint failed: F(10)\n";
        return false;
    }
    if (!close_to(solve_fast(100, window, radial_eps), 9.217221972, 5e-10)) {
        std::cerr << "Checkpoint failed: F(100)\n";
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

    if (options.run_checkpoints && !run_checkpoints(options.window, options.radial_eps)) {
        return 2;
    }

    const double ans = solve_fast(options.d, options.window, options.radial_eps);
    std::cout << std::fixed << std::setprecision(9) << ans << '\n';
    return 0;
}
