#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

struct KahanSum {
    long double sum = 0.0L;
    long double c = 0.0L;

    void add(long double value) {
        long double y = value - c;
        long double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
};

std::vector<long double> build_reciprocals(int n) {
    std::vector<long double> inv(n + 1, 0.0L);
    for (int i = 1; i <= n; ++i) {
        inv[i] = 1.0L / static_cast<long double>(i);
    }
    return inv;
}

std::vector<long double> build_g_values(int n, const std::vector<long double>& inv) {
    const int N = n - 1;
    std::vector<long double> G(N + 1, 0.0L);
    G[0] = 1.0L;
    for (int m = 1; m <= N; ++m) {
        long double term = 1.0L;
        long double sum = 1.0L;
        for (int c = 1; c <= m; ++c) {
            term *= static_cast<long double>(m - c + 1) * inv[N - c + 1];
            sum += term;
        }
        G[m] = sum;
    }
    return G;
}

long double compute_sum_sequential(int n,
                                   const std::vector<long double>& G,
                                   const std::vector<long double>& inv) {
    const int N = n - 1;
    KahanSum total;
    for (int a = 0; a <= N; ++a) {
        long double f = 1.0L;
        const int max_b = N - a;
        for (int b = 0; b <= max_b; ++b) {
            const int m = N - a - b;
            total.add(f * G[m]);
            if (b < max_b) {
                f *= static_cast<long double>(N - a - b) * inv[N - b];
            }
        }
    }
    return total.sum;
}

long double compute_sum_parallel(int n,
                                 const std::vector<long double>& G,
                                 const std::vector<long double>& inv,
                                 int threads) {
    const int N = n - 1;
    if (threads <= 1 || N < 200) {
        return compute_sum_sequential(n, G, inv);
    }
    const int max_threads = std::min(threads, N + 1);
    const int block = (N + 1 + max_threads - 1) / max_threads;
    std::vector<long double> partial(max_threads, 0.0L);
    std::vector<std::thread> pool;
    pool.reserve(max_threads);
    for (int t = 0; t < max_threads; ++t) {
        const int start = t * block;
        const int end = std::min(N + 1, start + block);
        if (start >= end) {
            continue;
        }
        pool.emplace_back([&, start, end, t]() {
            KahanSum local;
            for (int a = start; a < end; ++a) {
                long double f = 1.0L;
                const int max_b = N - a;
                for (int b = 0; b <= max_b; ++b) {
                    const int m = N - a - b;
                    local.add(f * G[m]);
                    if (b < max_b) {
                        f *= static_cast<long double>(N - a - b) * inv[N - b];
                    }
                }
            }
            partial[t] = local.sum;
        });
    }
    for (auto& th : pool) {
        th.join();
    }
    KahanSum total;
    for (long double value : partial) {
        total.add(value);
    }
    return total.sum;
}

long double compute_probability(int n, int threads) {
    if (n <= 1) {
        return 1.0L;
    }
    const int N = n - 1;
    const std::vector<long double> inv = build_reciprocals(N);
    const std::vector<long double> G = build_g_values(n, inv);
    const long double sum = compute_sum_parallel(n, G, inv, threads);
    return sum / (static_cast<long double>(n) * static_cast<long double>(n));
}

bool approx_equal(long double a, long double b, long double tol) {
    return std::fabs(a - b) <= tol;
}

bool run_validation() {
    bool ok = true;
    const long double p3 = compute_probability(3, 1);
    const long double p3_expected = 17.0L / 18.0L;
    if (!approx_equal(p3, p3_expected, 1e-12L)) {
        std::cerr << "Validation failed for n=3: got " << std::setprecision(15)
                  << p3 << " expected " << p3_expected
                  << "\n";
        ok = false;
    }
    const long double p10 = compute_probability(10, 1);
    const long double p10_expected = 0.6760292265L;
    if (!approx_equal(p10, p10_expected, 1e-9L)) {
        std::cerr << "Validation failed for n=10: got " << std::setprecision(15)
                  << p10 << " expected " << p10_expected
                  << "\n";
        ok = false;
    }
    if (ok) {
        std::cerr << "Validation checkpoints passed (n=3, n=10).\n";
    }
    return ok;
}

int main(int argc, char** argv) {
    int n = 20000;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) {
        threads = 1;
    }
    if (argc >= 2) {
        n = std::stoi(argv[1]);
    }
    if (argc >= 3) {
        threads = std::max(1, std::stoi(argv[2]));
    }

    const bool ok = run_validation();
    if (!ok) {
        return 1;
    }

    const long double result = compute_probability(n, threads);
    if (result < -1e-12L || result > 1.0L + 1e-12L) {
        std::cerr << "Probability out of range: " << std::setprecision(15)
                  << result << "\n";
    }
    std::cout << std::fixed << std::setprecision(10) << result << "\n";
    return 0;
}
