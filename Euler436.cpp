#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Options {
    bool run_checkpoints = true;
    unsigned requested_threads = 0U;
};

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    std::uint64_t parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<std::uint64_t>(c - '0');
        if (parsed > static_cast<std::uint64_t>(std::numeric_limits<unsigned>::max())) {
            return false;
        }
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
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

unsigned pick_thread_count(const unsigned requested) {
    if (requested > 0U) {
        return requested;
    }
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

long double cdf_last_draw_crossing_t(const long double t, const long double z) {
    if (z <= 0.0L) {
        return 0.0L;
    }
    if (z >= 1.0L) {
        return 1.0L;
    }
    if (t <= 0.0L) {
        return z;
    }
    if (t < z) {
        return 1.0L - (1.0L - z) * std::exp(t);
    }
    return std::exp(t) * (std::exp(-z) - 1.0L + z);
}

long double exact_probability() {
    // Closed form from integrating the exact joint density of (x, overshoot at 1).
    const long double e = std::exp(1.0L);
    return (1.0L + 14.0L * e - 5.0L * e * e) / 4.0L;
}

long double simpson_integral_cdf_segment(const long double a,
                                         const long double b,
                                         const long double z,
                                         const int n_even) {
    assert(n_even > 0 && (n_even % 2 == 0));
    if (b <= a) {
        return 0.0L;
    }
    const long double h = (b - a) / static_cast<long double>(n_even);
    long double sum = cdf_last_draw_crossing_t(a, z) + cdf_last_draw_crossing_t(b, z);
    for (int i = 1; i < n_even; ++i) {
        const long double v = a + h * static_cast<long double>(i);
        const long double w = (i & 1) ? 4.0L : 2.0L;
        sum += w * cdf_last_draw_crossing_t(v, z);
    }
    return sum * h / 3.0L;
}

long double simpson_integral_cdf(const long double t, const long double z, const int n_even) {
    assert(n_even > 0 && (n_even % 2 == 0));
    if (t <= 0.0L) {
        return 0.0L;
    }

    const long double split = std::min(t, std::max(0.0L, z));
    if (split <= 0.0L || split >= t) {
        return simpson_integral_cdf_segment(0.0L, t, z, n_even);
    }

    return simpson_integral_cdf_segment(0.0L, split, z, n_even) +
           simpson_integral_cdf_segment(split, t, z, n_even);
}

long double probability_integrand(const long double w, const long double x) {
    const long double t = 1.0L - w;
    const long double q_xw = std::exp(1.0L - x + w);  // Joint density of (x, overshoot at 1).
    const long double y_beats_x = 1.0L - cdf_last_draw_crossing_t(t, x);
    return q_xw * y_beats_x;
}

long double inner_integral_for_w(const long double w, const int n_even_x) {
    assert(n_even_x > 0 && (n_even_x % 2 == 0));
    if (w >= 1.0L) {
        return 0.0L;
    }

    const long double a = w;
    const long double b = 1.0L;
    const long double h = (b - a) / static_cast<long double>(n_even_x);
    long double sum = probability_integrand(w, a) + probability_integrand(w, b);
    for (int i = 1; i < n_even_x; ++i) {
        const long double x = a + h * static_cast<long double>(i);
        const long double wt = (i & 1) ? 4.0L : 2.0L;
        sum += wt * probability_integrand(w, x);
    }
    return sum * h / 3.0L;
}

long double numerical_probability_check(const int n_even_w,
                                        const int n_even_x,
                                        const unsigned thread_count) {
    assert(n_even_w > 0 && (n_even_w % 2 == 0));
    assert(n_even_x > 0 && (n_even_x % 2 == 0));

    std::vector<long double> inner(static_cast<std::size_t>(n_even_w) + 1ULL, 0.0L);
    const long double hw = 1.0L / static_cast<long double>(n_even_w);
    std::atomic<int> next_i(0);

    auto worker = [&]() {
        while (true) {
            const int i = next_i.fetch_add(1, std::memory_order_relaxed);
            if (i > n_even_w) {
                return;
            }
            const long double w = hw * static_cast<long double>(i);
            inner[static_cast<std::size_t>(i)] = inner_integral_for_w(w, n_even_x);
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(thread_count);
    for (unsigned i = 0U; i < thread_count; ++i) {
        pool.emplace_back(worker);
    }
    for (auto& th : pool) {
        th.join();
    }

    long double sum = inner.front() + inner.back();
    for (int i = 1; i < n_even_w; ++i) {
        const long double wt = (i & 1) ? 4.0L : 2.0L;
        sum += wt * inner[static_cast<std::size_t>(i)];
    }
    return sum * hw / 3.0L;
}

void run_checkpoints(const unsigned thread_count) {
    constexpr long double kEps = 1e-12L;

    // Basic boundary checks for the CDF.
    assert(std::fabs(cdf_last_draw_crossing_t(0.0L, 0.37L) - 0.37L) < kEps);
    assert(std::fabs(cdf_last_draw_crossing_t(0.42L, 0.0L) - 0.0L) < kEps);
    assert(std::fabs(cdf_last_draw_crossing_t(0.42L, 1.0L) - 1.0L) < kEps);

    // Renewal equation check: G(t,z) = ∫_0^t G(v,z)dv + max(0, z - t).
    const std::vector<long double> ts = {0.11L, 0.29L, 0.53L, 0.77L, 0.93L};
    const std::vector<long double> zs = {0.07L, 0.19L, 0.41L, 0.68L, 0.95L};
    for (const long double t : ts) {
        for (const long double z : zs) {
            const long double lhs = cdf_last_draw_crossing_t(t, z);
            const long double integral = simpson_integral_cdf(t, z, 2000);
            const long double rhs = integral + ((z > t) ? (z - t) : 0.0L);
            assert(std::fabs(lhs - rhs) < 5e-11L);
        }
    }

    // Independent numerical integration of the full 2D probability.
    const long double p_exact = exact_probability();
    const long double p_numeric = numerical_probability_check(600, 1200, thread_count);
    assert(std::fabs(p_numeric - p_exact) < 2e-9L);
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const unsigned thread_count = pick_thread_count(options.requested_threads);
    if (options.run_checkpoints) {
        run_checkpoints(thread_count);
    }

    const long double ans = exact_probability();
    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(ans) << '\n';
    return 0;
}
