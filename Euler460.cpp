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
    bool run_checkpoints = true;
    int dmax = 80;
    int max_band = 80;
    int band_step = 10;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
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
        if (parse_int_after_prefix(arg, "--d=", options.d)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--dmax=", options.dmax)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--max-band=", options.max_band)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--band-step=", options.band_step)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return options.d >= 1 && options.dmax >= 1 && options.max_band >= 10 && options.band_step >= 1;
}

double segment_time(const int dx, const int y0, const int y1) {
    if (y0 == y1) {
        return static_cast<double>(dx) / static_cast<double>(y0);
    }
    const double dy = static_cast<double>(y1 - y0);
    const double log_delta = std::fabs(std::log(static_cast<double>(y1)) - std::log(static_cast<double>(y0)));
    return std::hypot(static_cast<double>(dx), dy) * (log_delta / std::fabs(dy));
}

struct ProfileResult {
    double value = std::numeric_limits<double>::infinity();
};

ProfileResult solve_profiled(const int d, const int band, const int dmax) {
    const double radius = std::sqrt((d * 0.5) * (d * 0.5) + 1.0);

    std::vector<std::vector<int>> y_candidates(static_cast<std::size_t>(d) + 1U);
    int max_y = 1;
    for (int x = 0; x <= d; ++x) {
        double inside = radius * radius - (x - d * 0.5) * (x - d * 0.5);
        if (inside < 0.0) {
            inside = 0.0;
        }
        const double y_circle = std::sqrt(inside);
        const int yf = std::max(1, static_cast<int>(std::floor(y_circle)));
        const int yc = std::max(1, static_cast<int>(std::ceil(y_circle)));

        std::vector<int> ys;
        ys.reserve(static_cast<std::size_t>(2 * band + 4));
        ys.push_back(1);
        ys.push_back(yf);
        ys.push_back(yc);
        for (int b = 1; b <= band; ++b) {
            ys.push_back(std::max(1, yf - b));
            ys.push_back(std::max(1, yc + b));
        }

        std::sort(ys.begin(), ys.end());
        ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
        max_y = std::max(max_y, ys.back());
        y_candidates[static_cast<std::size_t>(x)] = std::move(ys);
    }

    std::vector<double> log_y(static_cast<std::size_t>(max_y) + 1U, 0.0);
    for (int y = 1; y <= max_y; ++y) {
        log_y[static_cast<std::size_t>(y)] = std::log(static_cast<double>(y));
    }

    const double inf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> dp(static_cast<std::size_t>(d) + 1U);
    for (int x = 0; x <= d; ++x) {
        dp[static_cast<std::size_t>(x)].assign(y_candidates[static_cast<std::size_t>(x)].size(), inf);
    }

    {
        const auto& ys0 = y_candidates[0];
        const auto it = std::lower_bound(ys0.begin(), ys0.end(), 1);
        dp[0][static_cast<std::size_t>(it - ys0.begin())] = 0.0;
    }

    const auto vertical_closure = [&](const int x) {
        std::vector<double>& arr = dp[static_cast<std::size_t>(x)];
        const std::vector<int>& ys = y_candidates[static_cast<std::size_t>(x)];
        const int n = static_cast<int>(ys.size());
        if (n <= 1) {
            return;
        }

        std::vector<double> out(static_cast<std::size_t>(n), inf);

        double best = inf;
        for (int i = 0; i < n; ++i) {
            best = std::min(best, arr[static_cast<std::size_t>(i)] -
                                      log_y[static_cast<std::size_t>(ys[static_cast<std::size_t>(i)])]);
            out[static_cast<std::size_t>(i)] =
                std::min(out[static_cast<std::size_t>(i)],
                         best + log_y[static_cast<std::size_t>(ys[static_cast<std::size_t>(i)])]);
        }

        best = inf;
        for (int i = n - 1; i >= 0; --i) {
            best = std::min(best, arr[static_cast<std::size_t>(i)] +
                                      log_y[static_cast<std::size_t>(ys[static_cast<std::size_t>(i)])]);
            out[static_cast<std::size_t>(i)] =
                std::min(out[static_cast<std::size_t>(i)],
                         best - log_y[static_cast<std::size_t>(ys[static_cast<std::size_t>(i)])]);
        }

        arr.swap(out);
    };

    vertical_closure(0);

    for (int x = 0; x <= d; ++x) {
        vertical_closure(x);

        const std::vector<int>& src_ys = y_candidates[static_cast<std::size_t>(x)];
        const std::vector<double>& src_dp = dp[static_cast<std::size_t>(x)];
        const int src_n = static_cast<int>(src_ys.size());

        for (int si = 0; si < src_n; ++si) {
            const double base = src_dp[static_cast<std::size_t>(si)];
            if (!std::isfinite(base)) {
                continue;
            }
            const int y0 = src_ys[static_cast<std::size_t>(si)];

            // Same-height horizontal motion can always be decomposed into unit-x moves.
            if (x + 1 <= d) {
                const std::vector<int>& ys_next = y_candidates[static_cast<std::size_t>(x + 1)];
                const auto it = std::lower_bound(ys_next.begin(), ys_next.end(), y0);
                if (it != ys_next.end() && *it == y0) {
                    const std::size_t ti = static_cast<std::size_t>(it - ys_next.begin());
                    dp[static_cast<std::size_t>(x + 1)][ti] =
                        std::min(dp[static_cast<std::size_t>(x + 1)][ti], base + 1.0 / static_cast<double>(y0));
                }
            }

            const int max_dx = std::min(dmax, d - x);
            for (int dx = 1; dx <= max_dx; ++dx) {
                const int nx = x + dx;
                const std::vector<int>& dst_ys = y_candidates[static_cast<std::size_t>(nx)];
                std::vector<double>& dst_dp = dp[static_cast<std::size_t>(nx)];
                const int dst_n = static_cast<int>(dst_ys.size());

                for (int ti = 0; ti < dst_n; ++ti) {
                    const int y1 = dst_ys[static_cast<std::size_t>(ti)];
                    if (y1 == y0 && dx > 1) {
                        continue;
                    }

                    const double cand = base + segment_time(dx, y0, y1);
                    if (cand < dst_dp[static_cast<std::size_t>(ti)]) {
                        dst_dp[static_cast<std::size_t>(ti)] = cand;
                    }
                }
            }
        }
    }

    vertical_closure(d);

    const std::vector<int>& ys_end = y_candidates[static_cast<std::size_t>(d)];
    const auto it_end = std::lower_bound(ys_end.begin(), ys_end.end(), 1);
    return ProfileResult{dp[static_cast<std::size_t>(d)][static_cast<std::size_t>(it_end - ys_end.begin())]};
}

double solve(const int d, const int dmax, const int max_band, const int band_step) {
    double previous = solve_profiled(d, band_step, dmax).value;
    bool stable_once = false;

    for (int band = 2 * band_step; band <= max_band; band += band_step) {
        const double current = solve_profiled(d, band, dmax).value;
        const double diff = std::fabs(current - previous);
        if (diff < 1e-12) {
            if (stable_once) {
                return current;
            }
            stable_once = true;
        } else {
            stable_once = false;
        }
        previous = current;
    }

    return previous;
}

bool close_to(const double a, const double b, const double eps) {
    return std::fabs(a - b) <= eps;
}

bool run_checkpoints(const Options& options) {
    const double f4 = solve(4, options.dmax, options.max_band, options.band_step);
    if (!close_to(f4, 2.960516287, 5e-10)) {
        std::cerr << "Checkpoint failed: F(4)\n";
        return false;
    }

    const double f10 = solve(10, options.dmax, options.max_band, options.band_step);
    if (!close_to(f10, 4.668187834, 5e-10)) {
        std::cerr << "Checkpoint failed: F(10)\n";
        return false;
    }

    const double f100 = solve(100, options.dmax, options.max_band, options.band_step);
    if (!close_to(f100, 9.217221972, 5e-10)) {
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

    if (options.run_checkpoints && !run_checkpoints(options)) {
        return 2;
    }

    const double ans = solve(options.d, options.dmax, options.max_band, options.band_step);
    std::cout << std::fixed << std::setprecision(9) << ans << '\n';
    return 0;
}
