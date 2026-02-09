#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

struct Options {
    int n = 10'000;
    int threads = 0;
    bool run_checkpoints = true;
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

u64 prefix_pairs(const int q, const int kmax) {
    const u64 qq = static_cast<u64>(q);
    const u64 km = static_cast<u64>(kmax + 1);
    return qq * km - (qq * (qq - 1ULL)) / 2ULL;
}

int find_boundary_for_target(const u64 target, const int kmax) {
    int lo = 0;
    int hi = kmax + 1;
    while (lo < hi) {
        const int mid = lo + (hi - lo) / 2;
        if (prefix_pairs(mid, kmax) >= target) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return lo;
}

std::vector<PairEntry> build_pairs(const std::vector<double>& values, int thread_count) {
    const int kmax = static_cast<int>(values.size()) - 1;
    const u64 total_pairs = prefix_pairs(kmax + 1, kmax);
    std::vector<PairEntry> pairs(static_cast<std::size_t>(total_pairs));

    if (thread_count <= 0) {
        thread_count = static_cast<int>(std::thread::hardware_concurrency());
    }
    if (thread_count <= 0) {
        thread_count = 4;
    }
    thread_count = std::min(thread_count, kmax + 1);

    std::vector<int> boundaries(static_cast<std::size_t>(thread_count + 1), 0);
    boundaries[0] = 0;
    boundaries[static_cast<std::size_t>(thread_count)] = kmax + 1;
    for (int t = 1; t < thread_count; ++t) {
        const u64 target =
            (total_pairs * static_cast<u64>(t)) / static_cast<u64>(thread_count);
        boundaries[static_cast<std::size_t>(t)] = find_boundary_for_target(target, kmax);
    }

    std::vector<u32> sq(static_cast<std::size_t>(kmax + 1), 0U);
    for (int i = 0; i <= kmax; ++i) {
        sq[static_cast<std::size_t>(i)] = static_cast<u32>(i) * static_cast<u32>(i);
    }

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));
    for (int t = 0; t < thread_count; ++t) {
        workers.emplace_back([&, t]() {
            const int ibegin = boundaries[static_cast<std::size_t>(t)];
            const int iend = boundaries[static_cast<std::size_t>(t + 1)];
            std::size_t out = static_cast<std::size_t>(prefix_pairs(ibegin, kmax));

            for (int i = ibegin; i < iend; ++i) {
                const double vi = values[static_cast<std::size_t>(i)];
                const u32 i_sq = sq[static_cast<std::size_t>(i)];
                for (int j = i; j <= kmax; ++j) {
                    pairs[out].sum = vi + values[static_cast<std::size_t>(j)];
                    pairs[out].sq = i_sq + sq[static_cast<std::size_t>(j)];
                    ++out;
                }
            }
        });
    }
    for (std::thread& w : workers) {
        w.join();
    }

    return pairs;
}

u64 solve(const int n, const int thread_count) {
    const long double pi = std::acos(-1.0L);

    const int kmax = compute_kmax(n);
    std::vector<double> values(static_cast<std::size_t>(kmax + 1), 0.0);
    for (int k = 0; k <= kmax; ++k) {
        values[static_cast<std::size_t>(k)] =
            static_cast<double>(std::expm1(static_cast<long double>(k) / n));
    }

    std::vector<PairEntry> pairs = build_pairs(values, thread_count);
    std::sort(pairs.begin(), pairs.end(), [](const PairEntry& a, const PairEntry& b) {
        if (a.sum < b.sum) {
            return true;
        }
        if (a.sum > b.sum) {
            return false;
        }
        return a.sq < b.sq;
    });

    std::size_t left = 0;
    std::size_t right = pairs.size() - 1;
    long double best_error = std::numeric_limits<long double>::infinity();
    u64 best_g = std::numeric_limits<u64>::max();

    while (left <= right) {
        const long double total =
            static_cast<long double>(pairs[left].sum) + static_cast<long double>(pairs[right].sum);
        const long double error = std::fabs(total - pi);
        const u64 g = static_cast<u64>(pairs[left].sq) + static_cast<u64>(pairs[right].sq);

        if (error + 1e-21L < best_error) {
            best_error = error;
            best_g = g;
        } else if (std::fabs(error - best_error) <= 1e-21L && g < best_g) {
            best_g = g;
        }

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
