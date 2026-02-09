#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 n = 10000000000ULL;
    int m = 4000;
    int threads = 0;
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
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = static_cast<u64>(std::stoull(tail));
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
        if (parse_u64_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--m=", options.m) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL && options.m >= 0 && options.threads >= 0;
}

int choose_thread_count(int requested, u64 work_items) {
    if (work_items <= 1ULL) {
        return 1;
    }
    int threads = requested;
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
    }
    if (threads <= 0) {
        threads = 4;
    }
    if (threads > static_cast<int>(work_items)) {
        threads = static_cast<int>(work_items);
    }
    if (threads < 1) {
        threads = 1;
    }
    return threads;
}

long double r_value(u64 n, u64 i) {
    // r_i = 2*P(not flipped in one turn) - 1.
    const long double nf = static_cast<long double>(n);
    const long double n2 = nf * nf;
    const long double a = static_cast<long double>(i - 1ULL);
    const long double b = static_cast<long double>(n - i);
    return (2.0L * (a * a + b * b) - n2) / n2;
}

long double expected_exact(u64 n, int m) {
    long double ans = 0.0L;
    for (u64 i = 1; i <= n; ++i) {
        const long double r = r_value(n, i);
        ans += 0.5L * (1.0L + std::powl(r, static_cast<long double>(m)));
    }
    return ans;
}

u64 edge_limit_by_rho(u64 n, long double rho) {
    u64 lo = 0;
    u64 hi = n / 2ULL;
    while (lo < hi) {
        const u64 mid = lo + (hi - lo + 1ULL) / 2ULL;
        if (r_value(n, mid) > rho) {
            lo = mid;
        } else {
            hi = mid - 1ULL;
        }
    }
    return lo;
}

long double edge_sum_parallel(u64 n, int m, u64 edge_count, int requested_threads) {
    if (edge_count == 0ULL) {
        return 0.0L;
    }
    const int thread_count = choose_thread_count(requested_threads, edge_count);
    std::vector<long double> partial(static_cast<std::size_t>(thread_count), 0.0L);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(thread_count));

    for (int t = 0; t < thread_count; ++t) {
        const u64 begin = edge_count * static_cast<u64>(t) / static_cast<u64>(thread_count) + 1ULL;
        const u64 end =
            edge_count * static_cast<u64>(t + 1) / static_cast<u64>(thread_count);
        workers.emplace_back([&, begin, end, t]() {
            long double local = 0.0L;
            for (u64 i = begin; i <= end; ++i) {
                const long double r = r_value(n, i);
                local += std::powl(r, static_cast<long double>(m));
            }
            partial[static_cast<std::size_t>(t)] = local;
        });
    }
    for (auto& w : workers) {
        w.join();
    }

    long double sum = 0.0L;
    for (long double x : partial) {
        sum += x;
    }
    return sum;
}

long double expected_fast(u64 n, int m, int threads) {
    if (m == 0) {
        return static_cast<long double>(n);
    }
    if (n <= 2000000ULL) {
        return expected_exact(n, m);
    }

    // Tail bound target: omitted contribution to E is < 1e-4.
    const long double eps = 1e-4L;
    const long double rho = std::expl(std::log((2.0L * eps) / static_cast<long double>(n)) /
                                      static_cast<long double>(m));

    const u64 l = edge_limit_by_rho(n, rho);
    const long double edge = edge_sum_parallel(n, m, l, threads);

    // E = N/2 + 1/2 * sum_i r_i^M.
    // Approximation keeps the 2*L edge terms exactly (by symmetry) and drops middle terms.
    long double estimate = static_cast<long double>(n) * 0.5L + edge;

    // Safety check for rounding to 2 decimals.
    const u64 omitted = n - 2ULL * l;
    const long double err_bound =
        0.5L * static_cast<long double>(omitted) * std::powl(rho, static_cast<long double>(m));
    if (err_bound >= 0.005L) {
        // Fallback: tighter threshold if needed.
        const long double rho2 = rho * 0.99L;
        const u64 l2 = edge_limit_by_rho(n, rho2);
        const long double edge2 = edge_sum_parallel(n, m, l2, threads);
        estimate = static_cast<long double>(n) * 0.5L + edge2;
    }

    return estimate;
}

bool close_to(long double a, long double b, long double tol) {
    return std::fabsl(a - b) <= tol;
}

bool run_checkpoints() {
    if (!close_to(expected_exact(3, 1), 10.0L / 9.0L, 1e-15L)) {
        std::cerr << "Checkpoint failed: E(3,1)\n";
        return false;
    }
    if (!close_to(expected_exact(3, 2), 5.0L / 3.0L, 1e-15L)) {
        std::cerr << "Checkpoint failed: E(3,2)\n";
        return false;
    }
    if (!close_to(expected_exact(10, 4), 5.15702608L, 1e-12L)) {
        std::cerr << "Checkpoint failed: E(10,4)\n";
        return false;
    }
    if (!close_to(expected_exact(100, 10), 51.8928010382861L, 1e-10L)) {
        std::cerr << "Checkpoint failed: E(100,10)\n";
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

    const long double answer = expected_fast(options.n, options.m, options.threads);
    std::cout << std::fixed << std::setprecision(2) << static_cast<double>(answer) << '\n';
    return 0;
}
