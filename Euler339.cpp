#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Options {
    int n = 10'000;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, std::uint64_t& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::uint64_t parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const std::uint64_t digit = static_cast<std::uint64_t>(c - '0');
        if (parsed > (std::numeric_limits<std::uint64_t>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    std::uint64_t parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<std::uint64_t>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        std::uint64_t parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            if (parsed_u64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
                std::cerr << "--n is too large for this build.\n";
                return false;
            }
            options.n = static_cast<int>(parsed_u64);
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n <= 0) {
        std::cerr << "--n must be >= 1.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

// Fast solver based on the optimal cap policy: after each bleat, if (white >= black),
// remove white sheep so that white = black - 1.
long double solve_fast(int n) {
    const int max_sum = 2 * n;

    std::vector<long double> prev(1ULL, 0.0L);
    std::vector<long double> cur;
    std::vector<long double> b;
    std::vector<long double> c;
    std::vector<long double> d;

    for (int s = 1; s <= max_sum; ++s) {
        cur.assign(static_cast<std::size_t>(s) + 1ULL, 0.0L);
        cur[0] = static_cast<long double>(s);
        cur[static_cast<std::size_t>(s)] = 0.0L;

        if (s >= 2) {
            const int t = (s - 1) / 2;  // continuation region: 1..t

            for (int w = t + 1; w <= s - 1; ++w) {
                cur[static_cast<std::size_t>(w)] = prev[static_cast<std::size_t>(w - 1)];
            }

            if (t >= 1) {
                const long double s_ld = static_cast<long double>(s);
                const long double right_boundary = prev[static_cast<std::size_t>(t)];

                b.assign(static_cast<std::size_t>(t) + 1ULL, 0.0L);
                c.assign(static_cast<std::size_t>(t) + 1ULL, 0.0L);
                d.assign(static_cast<std::size_t>(t) + 1ULL, 0.0L);

                b[1] = 1.0L;
                c[1] = (t >= 2) ? (-1.0L / s_ld) : 0.0L;
                d[1] = s_ld - 1.0L;
                if (t == 1) {
                    d[1] += right_boundary / s_ld;
                }

                for (int w = 2; w <= t; ++w) {
                    const long double a = -static_cast<long double>(s - w) / s_ld;
                    long double bw = 1.0L;
                    long double cw = (w < t) ? (-static_cast<long double>(w) / s_ld) : 0.0L;
                    long double dw = 0.0L;
                    if (w == t) {
                        dw += static_cast<long double>(w) * right_boundary / s_ld;
                    }

                    const long double mult = a / b[static_cast<std::size_t>(w - 1)];
                    bw -= mult * c[static_cast<std::size_t>(w - 1)];
                    dw -= mult * d[static_cast<std::size_t>(w - 1)];

                    b[static_cast<std::size_t>(w)] = bw;
                    c[static_cast<std::size_t>(w)] = cw;
                    d[static_cast<std::size_t>(w)] = dw;
                }

                cur[static_cast<std::size_t>(t)] =
                    d[static_cast<std::size_t>(t)] / b[static_cast<std::size_t>(t)];
                for (int w = t - 1; w >= 1; --w) {
                    cur[static_cast<std::size_t>(w)] =
                        (d[static_cast<std::size_t>(w)] -
                         c[static_cast<std::size_t>(w)] * cur[static_cast<std::size_t>(w + 1)]) /
                        b[static_cast<std::size_t>(w)];
                }
            }
        }

        prev.swap(cur);
    }

    const long double left = prev[static_cast<std::size_t>(n - 1)];
    const long double right = prev[static_cast<std::size_t>(n + 1)];
    return 0.5L * (left + right);
}

std::vector<std::vector<long double>> solve_exact_table(int max_sum,
                                                        long double tolerance = 1e-18L) {
    if (max_sum < 0) {
        throw std::runtime_error("max_sum must be non-negative.");
    }

    constexpr int kMaxIterations = 1'000'000;

    std::vector<std::vector<long double>> value(static_cast<std::size_t>(max_sum) + 1ULL);
    value[0] = {0.0L};

    for (int s = 1; s <= max_sum; ++s) {
        std::vector<long double> obstacle(static_cast<std::size_t>(s) + 1ULL, 0.0L);
        for (int w = 1; w < s; ++w) {
            obstacle[static_cast<std::size_t>(w)] =
                value[static_cast<std::size_t>(s - 1)][static_cast<std::size_t>(w - 1)];
        }

        std::vector<long double> cur(static_cast<std::size_t>(s) + 1ULL, 0.0L);
        cur[0] = static_cast<long double>(s);
        cur[static_cast<std::size_t>(s)] = 0.0L;
        for (int w = 1; w < s; ++w) {
            cur[static_cast<std::size_t>(w)] = obstacle[static_cast<std::size_t>(w)];
        }

        const long double s_ld = static_cast<long double>(s);
        for (int iter = 0; iter < kMaxIterations; ++iter) {
            long double max_delta = 0.0L;
            long double left_new = cur[0];
            for (int w = 1; w < s; ++w) {
                const std::size_t wi = static_cast<std::size_t>(w);
                const long double old_value = cur[wi];
                const long double continuation =
                    (static_cast<long double>(w) / s_ld) * cur[wi + 1ULL] +
                    (static_cast<long double>(s - w) / s_ld) * left_new;
                const long double new_value = std::max(obstacle[wi], continuation);
                cur[wi] = new_value;
                const long double delta = std::fabsl(new_value - old_value);
                max_delta = std::max(max_delta, delta);
                left_new = new_value;
            }

            if (max_delta < tolerance) {
                break;
            }
            if (iter + 1 == kMaxIterations) {
                throw std::runtime_error("Exact checkpoint solver did not converge.");
            }
        }

        value[static_cast<std::size_t>(s)] = std::move(cur);
    }

    return value;
}

long double expectation_from_table(const std::vector<std::vector<long double>>& table, int n) {
    const int s = 2 * n;
    const std::vector<long double>& diag = table[static_cast<std::size_t>(s)];
    return 0.5L * (diag[static_cast<std::size_t>(n - 1)] + diag[static_cast<std::size_t>(n + 1)]);
}

long double solve_exact_small(int n) {
    const auto table = solve_exact_table(2 * n);
    return expectation_from_table(table, n);
}

bool verify_cap_policy_on_small_domain(int max_sum, long double tolerance) {
    const auto table = solve_exact_table(max_sum);

    for (int b = 1; b <= max_sum; ++b) {
        for (int w = 0; w + b <= max_sum; ++w) {
            const int expected_best = std::min(w, b - 1);

            int best_k = 0;
            long double best_value = static_cast<long double>(b);  // k = 0

            for (int k = 1; k <= w; ++k) {
                const int s = k + b;
                const long double s_ld = static_cast<long double>(s);
                const long double candidate =
                    (static_cast<long double>(k) / s_ld) *
                        table[static_cast<std::size_t>(s)][static_cast<std::size_t>(k + 1)] +
                    (static_cast<long double>(b) / s_ld) *
                        table[static_cast<std::size_t>(s)][static_cast<std::size_t>(k - 1)];

                if (candidate > best_value + tolerance) {
                    best_value = candidate;
                    best_k = k;
                }
            }

            if (best_k != expected_best) {
                std::cerr << "Cap-policy checkpoint failed at state (w=" << w
                          << ", b=" << b << "). "
                          << "Expected best k=" << expected_best
                          << ", got k=" << best_k << '\n';
                return false;
            }
        }
    }

    return true;
}

bool run_checkpoints(const Options& options) {
    constexpr long double kGivenE5 = 6.871346L;
    constexpr long double kGivenTolerance = 0.5e-6L;
    constexpr long double kFastVsExactTolerance = 1e-11L;

    if (!verify_cap_policy_on_small_domain(80, 1e-12L)) {
        return false;
    }
    std::cout << "Checkpoint passed: cap policy verified for totals up to 80.\n";

    const std::vector<int> test_ns = {1, 5, 10, 20, 30};
    struct CheckResult {
        int n = 0;
        long double fast = 0.0L;
        long double exact = 0.0L;
    };

    std::vector<CheckResult> results(test_ns.size());
    const unsigned threads = choose_thread_count(
        options.allow_multithreading, options.requested_threads, test_ns.size());

    auto worker = [&](unsigned tid) {
        for (std::size_t i = tid; i < test_ns.size(); i += threads) {
            const int n = test_ns[i];
            results[i].n = n;
            results[i].fast = solve_fast(n);
            results[i].exact = solve_exact_small(n);
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threads > 0U ? threads - 1U : 0U);
    for (unsigned t = 1U; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    worker(0U);
    for (std::thread& thread : pool) {
        thread.join();
    }

    for (const CheckResult& result : results) {
        const long double delta = std::fabsl(result.fast - result.exact);
        if (delta > kFastVsExactTolerance) {
            std::cerr << std::setprecision(18);
            std::cerr << "Fast-vs-exact checkpoint failed for n=" << result.n
                      << ": fast=" << result.fast
                      << ", exact=" << result.exact
                      << ", |delta|=" << delta << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: n=" << result.n
                  << " fast/exact delta = " << std::setprecision(3) << std::scientific
                  << static_cast<double>(delta) << '\n';
        std::cout << std::defaultfloat;
    }

    const long double e5 = results[1].fast;  // n = 5
    if (std::fabsl(e5 - kGivenE5) > kGivenTolerance) {
        std::cerr << std::setprecision(12);
        std::cerr << "Given-value checkpoint failed: E(5) expected 6.871346 (rounded), got "
                  << e5 << '\n';
        return false;
    }
    std::cout << "Checkpoint passed: E(5) rounds to 6.871346.\n";

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    try {
        const auto start_time = std::chrono::steady_clock::now();

        if (options.run_checkpoints && !run_checkpoints(options)) {
            return 1;
        }

        const long double answer = solve_fast(options.n);

        const auto end_time = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end_time - start_time;

        std::cout << std::fixed << std::setprecision(6) << answer << '\n';
        std::cout << "Elapsed: " << std::setprecision(3) << std::defaultfloat
                  << elapsed.count() << " seconds\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
