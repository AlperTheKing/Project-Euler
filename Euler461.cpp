#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
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

int compute_kmax(const int n) {
    const long double pi = std::acos(-1.0L);
    int k = static_cast<int>(std::floor(static_cast<long double>(n) * std::log1pl(pi)));
    while (std::expm1((static_cast<long double>(k) + 1.0L) / static_cast<long double>(n)) <= pi) {
        ++k;
    }
    while (k > 0 && std::expm1(static_cast<long double>(k) / static_cast<long double>(n)) > pi) {
        --k;
    }
    return k;
}

inline void update_best(const long double error, const u64 g, long double& best_error, u64& best_g) {
    constexpr long double eps = 1e-21L;
    if (error + eps < best_error) {
        best_error = error;
        best_g = g;
    } else if (std::fabs(error - best_error) <= eps && g < best_g) {
        best_g = g;
    }
}

u64 solve(const int n, int thread_count) {
    const long double pi = std::acos(-1.0L);

    const int kmax = compute_kmax(n);
    std::vector<double> values(static_cast<std::size_t>(kmax + 1), 0.0);
    for (int k = 0; k <= kmax; ++k) {
        values[static_cast<std::size_t>(k)] =
            static_cast<double>(std::expm1(static_cast<long double>(k) / n));
    }

    std::vector<u32> sq(static_cast<std::size_t>(kmax + 1), 0U);
    for (int i = 0; i <= kmax; ++i) {
        sq[static_cast<std::size_t>(i)] = static_cast<u32>(i) * static_cast<u32>(i);
    }

    if (thread_count <= 0) {
        thread_count = 1;
    }
    (void)thread_count;

    constexpr double delta = 0.001;
    std::vector<PairEntry> small_sums;
    std::vector<PairEntry> large_sums;
    small_sums.reserve(20'000'000);
    large_sums.reserve(20'000'000);

    for (int a = 0; a <= kmax; ++a) {
        const double ea = values[static_cast<std::size_t>(a)];
        if (4.0 * ea > static_cast<double>(pi) + delta) {
            break;
        }
        const u32 a_sq = sq[static_cast<std::size_t>(a)];
        for (int b = a; b <= kmax; ++b) {
            const double eb = values[static_cast<std::size_t>(b)];
            if (ea + 3.0 * eb > static_cast<double>(pi) + delta) {
                break;
            }
            small_sums.push_back({ea + eb, a_sq + sq[static_cast<std::size_t>(b)]});
        }
    }

    for (int d = kmax; d >= 0; --d) {
        const double ed = values[static_cast<std::size_t>(d)];
        if (ed > static_cast<double>(pi) + delta) {
            continue;
        }
        if (4.0 * ed < static_cast<double>(pi) - delta) {
            break;
        }
        const u32 d_sq = sq[static_cast<std::size_t>(d)];
        for (int c = d; c >= 0; --c) {
            const double ec = values[static_cast<std::size_t>(c)];
            if (ec + ed > static_cast<double>(pi) + delta) {
                continue;
            }
            if (ed + 3.0 * ec < static_cast<double>(pi) - delta) {
                break;
            }
            large_sums.push_back({ec + ed, d_sq + sq[static_cast<std::size_t>(c)]});
        }
    }

    std::sort(small_sums.begin(), small_sums.end(), [](const PairEntry& a, const PairEntry& b) {
        return (a.sum < b.sum) || (a.sum == b.sum && a.sq < b.sq);
    });
    std::sort(large_sums.begin(), large_sums.end(), [](const PairEntry& a, const PairEntry& b) {
        return (a.sum < b.sum) || (a.sum == b.sum && a.sq < b.sq);
    });

    std::size_t left = 0;
    std::size_t right = large_sums.size() - 1;
    long double best_error = std::numeric_limits<long double>::infinity();
    u64 best_g = std::numeric_limits<u64>::max();

    while (left < small_sums.size()) {
        const long double total =
            static_cast<long double>(small_sums[left].sum) +
            static_cast<long double>(large_sums[right].sum);
        const long double error = std::fabs(total - pi);
        const u64 g = static_cast<u64>(small_sums[left].sq) + static_cast<u64>(large_sums[right].sq);
        update_best(error, g, best_error, best_g);

        if (total > pi) {
            if (right == 0U) {
                break;
            }
            --right;
        } else {
            ++left;
        }
    }
    return best_g;
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
