#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Options {
    int n = 10'000;
    int threads = 0;
    bool run_checkpoints = false;
};

struct PairEntry {
    double sum = 0.0;
    u32 sq = 0U;
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
        if (arg == "--checkpoints") {
            options.run_checkpoints = true;
            continue;
        }
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n <= 0) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    if (options.threads < 0) {
        std::cerr << "--threads must be non-negative.\n";
        return false;
    }
    return true;
}

std::optional<std::pair<int, int>> closest_sum(const std::vector<double>& x,
                                               const std::vector<double>& y,
                                               const double target,
                                               const double delta_start = 0.0001) {
    if (x.empty() || y.empty()) {
        return std::nullopt;
    }
    int best_i = -1;
    int best_j = -1;
    double best_delta = delta_start;

    int j = static_cast<int>(y.size()) - 1;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        if (j < 0) {
            break;
        }
        if (x[static_cast<std::size_t>(i)] + y[static_cast<std::size_t>(j)] < target - best_delta) {
            continue;
        }
        while (j >= 0 &&
               x[static_cast<std::size_t>(i)] + y[static_cast<std::size_t>(j)] > target &&
               x[static_cast<std::size_t>(i)] <= y[static_cast<std::size_t>(j)]) {
            --j;
        }
        if (j < 0 || x[static_cast<std::size_t>(i)] > y[static_cast<std::size_t>(j)]) {
            break;
        }

        int cand_j = j;
        double delta = target - x[static_cast<std::size_t>(i)] - y[static_cast<std::size_t>(j)];
        if (j < static_cast<int>(y.size()) - 1 &&
            x[static_cast<std::size_t>(i)] + y[static_cast<std::size_t>(j + 1)] > target &&
            x[static_cast<std::size_t>(i)] + y[static_cast<std::size_t>(j + 1)] - target < delta) {
            cand_j = j + 1;
            delta = x[static_cast<std::size_t>(i)] + y[static_cast<std::size_t>(cand_j)] - target;
        }

        if (delta < best_delta) {
            best_delta = delta;
            best_i = i;
            best_j = cand_j;
        }
    }
    if (best_i < 0 || best_j < 0) {
        return std::nullopt;
    }
    return std::make_pair(best_i, best_j);
}

u64 solve(const int n, int thread_count) {
    const double pi = std::acos(-1.0);
    const int bound = static_cast<int>(std::floor(static_cast<double>(n) * std::log(pi + 0.1)));
    std::vector<double> values(static_cast<std::size_t>(bound + 1), 0.0);
    for (int a = 0; a <= bound; ++a) {
        values[static_cast<std::size_t>(a)] = std::exp(static_cast<double>(a) / n) - 1.0;
    }

    if (thread_count <= 0) {
        thread_count = 1;
    }
    (void)thread_count;

    constexpr double delta = 0.001;
    std::vector<double> small_sums;
    std::vector<double> large_sums;
    small_sums.reserve(24'000'000);
    large_sums.reserve(24'000'000);

    for (int a = 0; a <= bound; ++a) {
        const double ea = values[static_cast<std::size_t>(a)];
        if (4.0 * ea > pi + delta) {
            break;
        }
        for (int b = a; b <= bound; ++b) {
            const double eb = values[static_cast<std::size_t>(b)];
            if (ea + 3.0 * eb > pi + delta) {
                break;
            }
            small_sums.push_back(ea + eb);
        }
    }

    for (int d = bound; d >= 0; --d) {
        const double ed = values[static_cast<std::size_t>(d)];
        if (ed > pi + delta) {
            continue;
        }
        if (4.0 * ed < pi - delta) {
            break;
        }
        for (int c = d; c >= 0; --c) {
            const double ec = values[static_cast<std::size_t>(c)];
            if (ec + ed > pi + delta) {
                continue;
            }
            if (ed + 3.0 * ec < pi - delta) {
                break;
            }
            large_sums.push_back(ec + ed);
        }
    }

    std::sort(small_sums.begin(), small_sums.end());
    std::sort(large_sums.begin(), large_sums.end());

    const auto ij_opt = closest_sum(small_sums, large_sums, pi);
    if (!ij_opt) {
        return 0ULL;
    }
    const int i = ij_opt->first;
    const int j = ij_opt->second;

    const auto ab_opt = closest_sum(values, values, small_sums[static_cast<std::size_t>(i)], 0.001);
    const auto cd_opt = closest_sum(values, values, large_sums[static_cast<std::size_t>(j)], 0.001);
    if (!ab_opt || !cd_opt) {
        return 0ULL;
    }

    std::vector<int> idx = {ab_opt->first, ab_opt->second, cd_opt->first, cd_opt->second};
    std::sort(idx.begin(), idx.end());
    return static_cast<u64>(idx[0]) * static_cast<u64>(idx[0]) +
           static_cast<u64>(idx[1]) * static_cast<u64>(idx[1]) +
           static_cast<u64>(idx[2]) * static_cast<u64>(idx[2]) +
           static_cast<u64>(idx[3]) * static_cast<u64>(idx[3]);
}

bool run_checkpoints(const int thread_count) {
    const u64 g200 = solve(200, thread_count);
    if (g200 != 64'658ULL) {
        std::cerr << "Checkpoint failed: g(200) expected 64658, got " << g200 << '\n';
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

    if (options.run_checkpoints && !run_checkpoints(options.threads)) {
        return 1;
    }

    std::cout << solve(options.n, options.threads) << '\n';
    return 0;
}
