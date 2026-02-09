#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 n = 1000000000000ULL;
    bool run_checkpoints = true;
};

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
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

u64 capacity_with_budget(const long double c, const long double a, const long double b, const u64 cap) {
    const int imax = static_cast<int>(std::floor(c / a));
    const int jmax = static_cast<int>(std::floor(c / b));
    const int cols = jmax + 2;

    std::vector<u64> dp(static_cast<std::size_t>((imax + 2) * cols), 0ULL);
    auto at = [&](const int i, const int j) -> u64& {
        return dp[static_cast<std::size_t>(i * cols + j)];
    };

    const long double eps = 1e-16L * (1.0L + c);
    for (int i = imax; i >= 0; --i) {
        const long double ai = static_cast<long double>(i) * a;
        for (int j = jmax; j >= 0; --j) {
            if (ai + static_cast<long double>(j) * b <= c + eps) {
                __uint128_t v = 1;
                v += at(i + 1, j);
                v += at(i, j + 1);
                if (v > cap) {
                    v = cap;
                }
                at(i, j) = static_cast<u64>(v);
            }
        }
    }

    return at(0, 0);
}

long double C_value(const u64 n, const long double a, const long double b) {
    long double lo = 0.0L;
    long double hi = 1.0L;
    while (capacity_with_budget(hi, a, b, n) < n) {
        hi *= 2.0L;
    }

    for (int it = 0; it < 90; ++it) {
        const long double mid = (lo + hi) / 2.0L;
        if (capacity_with_budget(mid, a, b, n) >= n) {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    return hi;
}

bool nearly_equal(const long double a, const long double b, const long double rel_tol) {
    const long double scale = std::max(std::fabsl(a), std::fabsl(b));
    if (scale == 0.0L) {
        return true;
    }
    return std::fabsl(a - b) <= rel_tol * scale;
}

bool run_checkpoints() {
    const long double c1 = C_value(5, 2.0L, 3.0L);
    if (!nearly_equal(c1, 5.0L, 1e-14L)) {
        std::cerr << "Checkpoint failed: C(5,2,3)\n";
        return false;
    }

    const long double c2 = C_value(500, std::sqrt(2.0L), std::sqrt(3.0L));
    if (!nearly_equal(c2, 13.22073197L, 5e-10L)) {
        std::cerr << "Checkpoint failed: C(500,sqrt2,sqrt3)\n";
        return false;
    }

    const long double c3 = C_value(20000, 5.0L, 7.0L);
    if (!nearly_equal(c3, 82.0L, 1e-14L)) {
        std::cerr << "Checkpoint failed: C(20000,5,7)\n";
        return false;
    }

    const long double c4 = C_value(2000000, std::sqrt(5.0L), std::sqrt(7.0L));
    if (!nearly_equal(c4, 49.63755955L, 5e-10L)) {
        std::cerr << "Checkpoint failed: C(2000000,sqrt5,sqrt7)\n";
        return false;
    }

    return true;
}

long double solve_sum(const u64 n) {
    std::vector<u64> fib(31, 0ULL);
    fib[1] = 1ULL;
    fib[2] = 1ULL;
    for (int k = 3; k <= 30; ++k) {
        fib[static_cast<std::size_t>(k)] = fib[static_cast<std::size_t>(k - 1)] + fib[static_cast<std::size_t>(k - 2)];
    }

    long double ans = 0.0L;
    for (int k = 1; k <= 30; ++k) {
        const long double a = std::sqrt(static_cast<long double>(k));
        const long double b = std::sqrt(static_cast<long double>(fib[static_cast<std::size_t>(k)]));
        ans += C_value(n, a, b);
    }
    return ans;
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

    const long double answer = solve_sum(options.n);
    std::cout << std::fixed << std::setprecision(8) << static_cast<double>(answer) << '\n';
    return 0;
}
