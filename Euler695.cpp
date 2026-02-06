#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;

constexpr int kDefaultBaseN = 8192;

struct Options {
    int base_n = kDefaultBaseN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

inline double median3(const double a, const double b, const double c) {
    return a + b + c - std::min(a, std::min(b, c)) - std::max(a, std::max(b, c));
}

// M(x, y): expectation over the 6 relative label permutations of the median area,
// after normalizing gaps so that a+b=1 and c+d=1.
inline double mean_median_normalized(const double x, const double y) {
    const double x0 = 1.0 - x;
    const double y0 = 1.0 - y;

    const double m1 = median3(x * y, x0 * y0, 1.0);
    const double m2 = median3(x * y, x0, y0);
    const double m3 = median3(x * y0, x0 * y, 1.0);
    const double m4 = median3(x * y0, x0, y);
    const double m5 = median3(x, x0 * y, y0);
    const double m6 = median3(x, x0 * y0, y);

    return (m1 + m2 + m3 + m4 + m5 + m6) * (1.0 / 6.0);
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
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
        parsed = parsed * 10ULL + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
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

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--n=", parsed_u32)) {
            if (parsed_u32 > static_cast<u32>(std::numeric_limits<int>::max())) {
                std::cerr << "--n is too large.\n";
                return false;
            }
            options.base_n = static_cast<int>(parsed_u32);
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

    if (options.base_n < 2) {
        std::cerr << "--n must be at least 2.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const int workload) {
    if (!allow_multithreading || workload < 2) {
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

struct NodesWeights {
    std::vector<double> x;
    std::vector<double> w;
};

NodesWeights gauss_legendre_unit_interval(const int n) {
    NodesWeights nw;
    nw.x.assign(static_cast<std::size_t>(n), 0.0);
    nw.w.assign(static_cast<std::size_t>(n), 0.0);

    const double pi = std::acos(-1.0);
    const double eps = 1e-15;
    const int half = (n + 1) / 2;

    for (int i = 0; i < half; ++i) {
        double z = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
        double z_prev = 0.0;
        double p1 = 0.0;
        double p2 = 0.0;
        double pp = 0.0;

        do {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * static_cast<double>(j) - 1.0) * z * p2 -
                      (static_cast<double>(j) - 1.0) * p3) /
                     static_cast<double>(j);
            }
            pp = static_cast<double>(n) * (z * p1 - p2) / (z * z - 1.0);
            z_prev = z;
            z = z_prev - p1 / pp;
        } while (std::abs(z - z_prev) > eps);

        const double weight = 2.0 / ((1.0 - z * z) * pp * pp);
        const int i_left = i;
        const int i_right = n - 1 - i;

        nw.x[static_cast<std::size_t>(i_left)] = 0.5 * (-z + 1.0);
        nw.x[static_cast<std::size_t>(i_right)] = 0.5 * (z + 1.0);
        nw.w[static_cast<std::size_t>(i_left)] = 0.5 * weight;
        nw.w[static_cast<std::size_t>(i_right)] = 0.5 * weight;
    }

    return nw;
}

double integrate_j(const int n,
                   const bool allow_multithreading,
                   const unsigned requested_threads) {
    const NodesWeights nw = gauss_legendre_unit_interval(n);

    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, n);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    std::vector<double> partial(threads, 0.0);
    std::atomic<int> next_i(0);

    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            double local_sum = 0.0;
            while (true) {
                const int i = next_i.fetch_add(1, std::memory_order_relaxed);
                if (i >= n) {
                    break;
                }

                const double xi = nw.x[static_cast<std::size_t>(i)];
                const double wi = nw.w[static_cast<std::size_t>(i)];

                double row = wi * wi * mean_median_normalized(xi, xi);
                for (int j = i + 1; j < n; ++j) {
                    row += 2.0 * wi * nw.w[static_cast<std::size_t>(j)] *
                           mean_median_normalized(xi, nw.x[static_cast<std::size_t>(j)]);
                }

                local_sum += row;
            }
            partial[static_cast<std::size_t>(t)] = local_sum;
        });
    }

    for (auto& th : pool) {
        th.join();
    }

    double total = 0.0;
    for (const double s : partial) {
        total += s;
    }
    return total;
}

bool nearly_equal(const double a, const double b, const double eps) {
    return std::abs(a - b) <= eps;
}

bool run_checkpoints() {
    {
        const double c1 = mean_median_normalized(0.0, 0.0);
        const double c2 = mean_median_normalized(0.0, 0.5);
        const double c3 = mean_median_normalized(0.25, 0.25);
        const double c4 = mean_median_normalized(0.0, 0.25);

        if (!nearly_equal(c1, 1.0 / 3.0, 1e-15)) {
            std::cerr << "Checkpoint failed: M(0,0) != 1/3.\n";
            return false;
        }
        if (!nearly_equal(c2, 0.5, 1e-15)) {
            std::cerr << "Checkpoint failed: M(0,0.5) != 1/2.\n";
            return false;
        }
        if (!nearly_equal(c3, 0.375, 1e-15)) {
            std::cerr << "Checkpoint failed: M(0.25,0.25) != 3/8.\n";
            return false;
        }
        if (!nearly_equal(c4, 5.0 / 12.0, 1e-15)) {
            std::cerr << "Checkpoint failed: M(0,0.25) != 5/12.\n";
            return false;
        }
    }

    {
        constexpr std::array<std::pair<double, double>, 4> points = {
            std::pair<double, double>{0.19, 0.73},
            std::pair<double, double>{0.41, 0.22},
            std::pair<double, double>{0.8, 0.3},
            std::pair<double, double>{0.63, 0.63},
        };

        for (const auto& [x, y] : points) {
            const double f = mean_median_normalized(x, y);
            const double fx = mean_median_normalized(1.0 - x, y);
            const double fy = mean_median_normalized(x, 1.0 - y);
            const double fs = mean_median_normalized(y, x);
            if (!nearly_equal(f, fx, 2e-15) || !nearly_equal(f, fy, 2e-15) ||
                !nearly_equal(f, fs, 2e-15)) {
                std::cerr << "Checkpoint failed: symmetry mismatch at x=" << x << ", y=" << y
                          << ".\n";
                return false;
            }
        }
    }

    {
        const double e256 = integrate_j(256, false, 1U) / 4.0;
        const double e512 = integrate_j(512, false, 1U) / 4.0;
        const double e1024 = integrate_j(1024, false, 1U) / 4.0;

        if (!nearly_equal(e256, 0.10177862574564227, 5e-14)) {
            std::cerr << "Checkpoint failed: n=256 quadrature mismatch.\n";
            return false;
        }
        if (!nearly_equal(e512, 0.10177867340598332, 5e-14)) {
            std::cerr << "Checkpoint failed: n=512 quadrature mismatch.\n";
            return false;
        }
        if (!nearly_equal(e1024, 0.10177868278503224, 5e-14)) {
            std::cerr << "Checkpoint failed: n=1024 quadrature mismatch.\n";
            return false;
        }
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
        return 1;
    }

    const int base_n = options.base_n;
    const int fine_n = base_n * 2;

    const auto t0 = std::chrono::steady_clock::now();

    const double e_base = integrate_j(base_n, options.allow_multithreading, options.requested_threads) / 4.0;
    const double e_fine = integrate_j(fine_n, options.allow_multithreading, options.requested_threads) / 4.0;

    // Leading quadrature error behaves like O(n^-2) because of piecewise-smooth boundaries.
    const double e_richardson = e_fine + (e_fine - e_base) / 3.0;

    const auto t1 = std::chrono::steady_clock::now();
    const double elapsed_sec = std::chrono::duration<double>(t1 - t0).count();

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "E_base(n=" << base_n << ") = " << e_base << '\n';
    std::cout << "E_fine(n=" << fine_n << ") = " << e_fine << '\n';
    std::cout << "E_richardson = " << e_richardson << '\n';
    std::cout << "Answer (10 d.p.) = " << std::setprecision(10) << e_richardson << '\n';
    std::cout << std::setprecision(6) << "Elapsed seconds = " << elapsed_sec << '\n';

    return 0;
}
