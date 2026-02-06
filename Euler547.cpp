#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr int kDefaultN = 40;
constexpr int kDefaultQuadratureOrder = 64;

struct Options {
    int n = kDefaultN;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
    int quadrature_order = kDefaultQuadratureOrder;
};

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    std::int64_t parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const int digit = c - '0';
        if (parsed > (std::numeric_limits<std::int64_t>::max() - digit) / 10) {
            return false;
        }
        parsed = parsed * 10 + digit;
    }

    if (parsed > std::numeric_limits<int>::max()) {
        return false;
    }

    value = static_cast<int>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    int parsed = 0;
    if (!parse_int_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed < 0) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        int n = 0;
        if (parse_int_after_prefix(arg, "--n=", n)) {
            options.n = n;
            continue;
        }

        int q = 0;
        if (parse_int_after_prefix(arg, "--quad=", q)) {
            options.quadrature_order = q;
            continue;
        }

        unsigned t = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", t)) {
            options.requested_threads = t;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.n < 3) {
        std::cerr << "n must be at least 3.\n";
        return false;
    }
    if (options.quadrature_order < 8 || (options.quadrature_order % 2) != 0) {
        std::cerr << "--quad must be an even integer >= 8.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    threads = std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
    return threads;
}

void legendre_polynomial_and_derivative(int n, long double x, long double& pn, long double& dpn) {
    long double pnm1 = 1.0L;
    long double pn_local = x;
    if (n == 0) {
        pn = 1.0L;
        dpn = 0.0L;
        return;
    }
    if (n == 1) {
        pn = x;
        dpn = 1.0L;
        return;
    }

    for (int k = 2; k <= n; ++k) {
        const long double pk =
            ((2.0L * k - 1.0L) * x * pn_local - (k - 1.0L) * pnm1) / static_cast<long double>(k);
        pnm1 = pn_local;
        pn_local = pk;
    }

    pn = pn_local;
    dpn = static_cast<long double>(n) * (x * pn_local - pnm1) / (x * x - 1.0L);
}

void build_gauss_legendre_01(int order,
                             std::vector<long double>& nodes,
                             std::vector<long double>& weights) {
    nodes.assign(order, 0.0L);
    weights.assign(order, 0.0L);

    const int half = (order + 1) / 2;
    const long double pi = std::acos(-1.0L);

    for (int i = 0; i < half; ++i) {
        long double z =
            std::cos(pi * (static_cast<long double>(i) + 0.75L) / (static_cast<long double>(order) + 0.5L));
        for (int iter = 0; iter < 80; ++iter) {
            long double pn = 0.0L;
            long double dpn = 0.0L;
            legendre_polynomial_and_derivative(order, z, pn, dpn);
            const long double dz = pn / dpn;
            z -= dz;
            if (std::abs(dz) < 1e-25L) {
                break;
            }
        }

        long double pn = 0.0L;
        long double dpn = 0.0L;
        legendre_polynomial_and_derivative(order, z, pn, dpn);
        const long double w = 2.0L / ((1.0L - z * z) * dpn * dpn);

        const int j = order - 1 - i;
        nodes[i] = (1.0L - z) * 0.5L;
        nodes[j] = (1.0L + z) * 0.5L;
        weights[i] = w * 0.5L;
        weights[j] = w * 0.5L;
    }
}

std::vector<std::vector<long double>> build_kernel_table(int max_n, int quadrature_order) {
    const int max_diff = max_n - 1;
    std::vector<std::vector<long double>> kernel(max_diff + 1,
                                                 std::vector<long double>(max_diff + 1, 0.0L));

    std::vector<long double> nodes;
    std::vector<long double> weights;
    build_gauss_legendre_01(quadrature_order, nodes, weights);

    for (int dx = 0; dx <= max_diff; ++dx) {
        for (int dy = 0; dy <= max_diff; ++dy) {
            long double sum = 0.0L;
            for (int i = 0; i < quadrature_order; ++i) {
                const long double u = nodes[i];
                const long double wu = weights[i];
                const long double ux = 1.0L - u;
                for (int j = 0; j < quadrature_order; ++j) {
                    const long double v = nodes[j];
                    const long double wv = weights[j];
                    const long double vy = 1.0L - v;

                    const long double p1 = std::hypotl(static_cast<long double>(dx) + u,
                                                       static_cast<long double>(dy) + v);
                    const long double p2 = std::hypotl(static_cast<long double>(dx) + u,
                                                       static_cast<long double>(dy) - v);
                    const long double p3 = std::hypotl(static_cast<long double>(dx) - u,
                                                       static_cast<long double>(dy) + v);
                    const long double p4 = std::hypotl(static_cast<long double>(dx) - u,
                                                       static_cast<long double>(dy) - v);
                    sum += wu * wv * ux * vy * (p1 + p2 + p3 + p4);
                }
            }
            kernel[dx][dy] = sum;
        }
    }

    return kernel;
}

inline long double kernel_at(const std::vector<std::vector<long double>>& kernel, int dx, int dy) {
    return kernel[std::abs(dx)][std::abs(dy)];
}

std::vector<std::vector<long double>> build_rect_self_table(
    int max_n,
    const std::vector<std::vector<long double>>& kernel) {
    std::vector<std::vector<long double>> self(max_n + 1,
                                               std::vector<long double>(max_n + 1, 0.0L));
    for (int w = 1; w <= max_n; ++w) {
        for (int h = 1; h <= max_n; ++h) {
            long double sum = 0.0L;
            for (int dx = -(w - 1); dx <= (w - 1); ++dx) {
                const long double cx = static_cast<long double>(w - std::abs(dx));
                for (int dy = -(h - 1); dy <= (h - 1); ++dy) {
                    const long double cy = static_cast<long double>(h - std::abs(dy));
                    sum += cx * cy * kernel_at(kernel, dx, dy);
                }
            }
            self[w][h] = sum;
        }
    }
    return self;
}

std::vector<std::vector<long double>> build_outer_prefix(
    int n,
    const std::vector<std::vector<long double>>& kernel) {
    std::vector<std::vector<long double>> w(n, std::vector<long double>(n, 0.0L));
    for (int hx = 0; hx < n; ++hx) {
        for (int hy = 0; hy < n; ++hy) {
            long double sum = 0.0L;
            for (int ox = 0; ox < n; ++ox) {
                for (int oy = 0; oy < n; ++oy) {
                    sum += kernel_at(kernel, ox - hx, oy - hy);
                }
            }
            w[hx][hy] = sum;
        }
    }

    std::vector<std::vector<long double>> prefix(n + 1, std::vector<long double>(n + 1, 0.0L));
    for (int x = 0; x < n; ++x) {
        long double row_acc = 0.0L;
        for (int y = 0; y < n; ++y) {
            row_acc += w[x][y];
            prefix[x + 1][y + 1] = prefix[x][y + 1] + row_acc;
        }
    }
    return prefix;
}

inline long double rect_sum(const std::vector<std::vector<long double>>& prefix,
                            int x0,
                            int y0,
                            int w,
                            int h) {
    const int x1 = x0 + w;
    const int y1 = y0 + h;
    return prefix[x1][y1] - prefix[x0][y1] - prefix[x1][y0] + prefix[x0][y0];
}

long double compute_s_for_n(int n,
                            const std::vector<std::vector<long double>>& kernel,
                            const std::vector<std::vector<long double>>& rect_self,
                            bool allow_multithreading,
                            unsigned requested_threads) {
    if (n < 3) {
        return 0.0L;
    }

    const std::vector<std::vector<long double>> prefix = build_outer_prefix(n, kernel);
    const long double f_outer_outer = rect_self[n][n];

    std::vector<std::pair<int, int>> tasks;
    tasks.reserve(static_cast<std::size_t>(n - 2) * static_cast<std::size_t>(n - 2));
    for (int x = 1; x <= n - 2; ++x) {
        for (int y = 1; y <= n - 2; ++y) {
            tasks.emplace_back(x, y);
        }
    }

    const unsigned thread_count =
        choose_thread_count(allow_multithreading, requested_threads, tasks.size());
    std::vector<long double> partial(thread_count, 0.0L);
    std::atomic<std::size_t> next_task(0);

    auto worker = [&](unsigned tid) {
        long double local_sum = 0.0L;
        while (true) {
            const std::size_t idx = next_task.fetch_add(1, std::memory_order_relaxed);
            if (idx >= tasks.size()) {
                break;
            }

            const int x = tasks[idx].first;
            const int y = tasks[idx].second;
            const long double f_hole_hole = rect_self[x][y];
            const long double area_ring =
                static_cast<long double>(n) * static_cast<long double>(n) -
                static_cast<long double>(x) * static_cast<long double>(y);
            const long double inv_area_sq = 1.0L / (area_ring * area_ring);

            for (int a = 1; a <= n - x - 1; ++a) {
                for (int b = 1; b <= n - y - 1; ++b) {
                    const long double f_outer_hole = rect_sum(prefix, a, b, x, y);
                    const long double f_ring_ring =
                        f_outer_outer - 2.0L * f_outer_hole + f_hole_hole;
                    local_sum += f_ring_ring * inv_area_sq;
                }
            }
        }
        partial[tid] = local_sum;
    };

    std::vector<std::thread> threads;
    threads.reserve(thread_count > 0 ? thread_count - 1 : 0);
    for (unsigned t = 1; t < thread_count; ++t) {
        threads.emplace_back(worker, t);
    }
    worker(0);
    for (auto& th : threads) {
        th.join();
    }

    long double total = 0.0L;
    for (long double v : partial) {
        total += v;
    }
    return total;
}

long long rounded4(long double x) { return std::llround(x * 10000.0L); }

bool run_checkpoints(const std::vector<std::vector<long double>>& kernel,
                     const std::vector<std::vector<long double>>& rect_self) {
    // Checkpoint 1: unit-square mean distance.
    const long double sqrt2 = std::sqrt(2.0L);
    const long double expected_k00 =
        (2.0L + sqrt2 + 5.0L * std::log(1.0L + sqrt2)) / 15.0L;
    const long double got_k00 = kernel[0][0];
    if (std::abs(got_k00 - expected_k00) > 1e-10L) {
        std::cerr << "Checkpoint failed: kernel(0,0) mismatch. got=" << std::setprecision(16)
                  << got_k00 << " expected=" << expected_k00 << '\n';
        return false;
    }

    // Checkpoint 2: n=4 has 9 hollow laminae.
    long long laminae_count = 0;
    for (int x = 1; x <= 2; ++x) {
        for (int y = 1; y <= 2; ++y) {
            laminae_count += static_cast<long long>(4 - x - 1) * static_cast<long long>(4 - y - 1);
        }
    }
    if (laminae_count != 9LL) {
        std::cerr << "Checkpoint failed: laminae count for n=4 is " << laminae_count
                  << ", expected 9\n";
        return false;
    }

    // Checkpoint 3: statement samples.
    const long double s3 = compute_s_for_n(3, kernel, rect_self, false, 1);
    const long double s4 = compute_s_for_n(4, kernel, rect_self, false, 1);
    if (rounded4(s3) != 16514LL) {
        std::cerr << "Checkpoint failed: S(3) rounded to 4 decimals is " << std::fixed
                  << std::setprecision(4) << s3 << ", expected 1.6514\n";
        return false;
    }
    if (rounded4(s4) != 196564LL) {
        std::cerr << "Checkpoint failed: S(4) rounded to 4 decimals is " << std::fixed
                  << std::setprecision(4) << s4 << ", expected 19.6564\n";
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

    const int max_n = std::max(options.n, 4);
    const std::vector<std::vector<long double>> kernel =
        build_kernel_table(max_n, options.quadrature_order);
    const std::vector<std::vector<long double>> rect_self = build_rect_self_table(max_n, kernel);

    if (options.run_checkpoints && !run_checkpoints(kernel, rect_self)) {
        return 1;
    }

    const long double answer = compute_s_for_n(
        options.n, kernel, rect_self, options.allow_multithreading, options.requested_threads);
    std::cout << std::fixed << std::setprecision(4) << answer << '\n';
    return 0;
}

