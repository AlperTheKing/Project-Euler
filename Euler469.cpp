#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;

struct Options {
    u64 n = 1'000'000'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
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
    if (options.n == 0ULL) {
        std::cerr << "--n must be positive.\n";
        return false;
    }
    return true;
}

long double alt_exp_partial(const u64 m) {
    long double term = 1.0L;
    long double sum = 1.0L;
    for (u64 k = 1; k <= m; ++k) {
        term *= -2.0L / static_cast<long double>(k);
        sum += term;
    }
    return sum;
}

long double exact_cycle_expectation_small(const u64 n) {
    if (n == 1ULL) {
        return 0.0L;
    }
    if (n == 2ULL) {
        return 0.5L;
    }
    if (n == 3ULL) {
        return 2.0L / 3.0L;
    }

    const u64 m = n - 3ULL;
    const long double e_m = alt_exp_partial(m);
    const long double e_m1 = alt_exp_partial(m - 1ULL);

    // For path length m: A(m) = (m+1)/2 - ((m+1)/2)E_m - E_{m-1},
    // where E_t = sum_{k=0}^t (-2)^k/k!.
    const long double a_m =
        (static_cast<long double>(m + 1ULL) * 0.5L) -
        (static_cast<long double>(m + 1ULL) * 0.5L) * e_m - e_m1;

    // For cycle length n: occupied expectation B(n)=1+A(n-3), so empty fraction is 1-B(n)/n.
    return 1.0L - (1.0L + a_m) / static_cast<long double>(n);
}

long double solve(const u64 n) {
    // For large n, E(n) = (1 + exp(-2))/2 up to a remainder from a truncated
    // alternating exponential tail. At n=1e18 this error is far below 1e-14.
    if (n > 200ULL) {
        return 0.5L * (1.0L + std::expl(-2.0L));
    }
    return exact_cycle_expectation_small(n);
}

bool almost_equal(const long double a, const long double b, const long double eps = 1e-18L) {
    return std::fabsl(a - b) <= eps;
}

bool run_checkpoints() {
    if (!almost_equal(solve(4ULL), 0.5L)) {
        std::cerr << "Checkpoint failed: E(4)\n";
        return false;
    }
    if (!almost_equal(solve(6ULL), 5.0L / 9.0L)) {
        std::cerr << "Checkpoint failed: E(6)\n";
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
        return 1;
    }

    std::cout << std::fixed << std::setprecision(14) << solve(options.n) << '\n';
    return 0;
}
