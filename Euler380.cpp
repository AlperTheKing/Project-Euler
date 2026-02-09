#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Options {
    int m = 100;
    int n = 500;
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

    int parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--m=", options.m) ||
            parse_int_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.m >= 1 && options.n >= 1;
}

long double compute_log10_spanning_trees(const int m, const int n) {
    const long double pi = acosl(-1.0L);

    std::vector<long double> cos_m(static_cast<std::size_t>(m));
    std::vector<long double> cos_n(static_cast<std::size_t>(n));
    for (int i = 0; i < m; ++i) {
        cos_m[static_cast<std::size_t>(i)] = cosl(pi * static_cast<long double>(i) / m);
    }
    for (int j = 0; j < n; ++j) {
        cos_n[static_cast<std::size_t>(j)] = cosl(pi * static_cast<long double>(j) / n);
    }

    // Kahan summation keeps rounding error tiny across many logarithms.
    long double sum_log10 = 0.0L;
    long double compensation = 0.0L;
    for (int i = 0; i < m; ++i) {
        const long double ci = cos_m[static_cast<std::size_t>(i)];
        for (int j = 0; j < n; ++j) {
            if (i == 0 && j == 0) {
                continue;
            }
            const long double cj = cos_n[static_cast<std::size_t>(j)];
            const long double eigenvalue = 4.0L - 2.0L * ci - 2.0L * cj;
            const long double y = log10l(eigenvalue) - compensation;
            const long double t = sum_log10 + y;
            compensation = (t - sum_log10) - y;
            sum_log10 = t;
        }
    }

    sum_log10 -= log10l(static_cast<long double>(m) * static_cast<long double>(n));
    return sum_log10;
}

std::string format_scientific_from_log10(const long double log10_value, const int significant_digits) {
    const long double exponent_ld = floorl(log10_value);
    long long exponent = static_cast<long long>(exponent_ld);
    long double mantissa = powl(10.0L, log10_value - exponent_ld);

    long double scale = 1.0L;
    for (int i = 1; i < significant_digits; ++i) {
        scale *= 10.0L;
    }

    long double rounded = roundl(mantissa * scale);
    if (rounded >= 10.0L * scale) {
        rounded /= 10.0L;
        ++exponent;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(significant_digits - 1)
        << static_cast<double>(rounded / scale) << 'e' << exponent;
    return oss.str();
}

bool is_nearly_equal(const long double a, const long double b, const long double rel_tol) {
    const long double diff = fabsl(a - b);
    const long double scale = std::max(fabsl(a), fabsl(b));
    if (scale == 0.0L) {
        return diff <= rel_tol;
    }
    return diff <= rel_tol * scale;
}

bool run_checkpoints() {
    if (!is_nearly_equal(compute_log10_spanning_trees(1, 1), 0.0L, 1e-18L)) {
        std::cerr << "Checkpoint failed: C(1,1) != 1\n";
        return false;
    }

    const long double c22 = powl(10.0L, compute_log10_spanning_trees(2, 2));
    if (llround(c22) != 4LL) {
        std::cerr << "Checkpoint failed: C(2,2)\n";
        return false;
    }

    const long double c34 = powl(10.0L, compute_log10_spanning_trees(3, 4));
    if (llround(c34) != 2415LL) {
        std::cerr << "Checkpoint failed: C(3,4)\n";
        return false;
    }

    const std::string c912 = format_scientific_from_log10(compute_log10_spanning_trees(9, 12), 5);
    if (c912 != "2.5720e46") {
        std::cerr << "Checkpoint failed: C(9,12), got " << c912 << '\n';
        return false;
    }

    const long double symmetry_lhs = compute_log10_spanning_trees(3, 4);
    const long double symmetry_rhs = compute_log10_spanning_trees(4, 3);
    if (!is_nearly_equal(symmetry_lhs, symmetry_rhs, 1e-16L)) {
        std::cerr << "Checkpoint failed: symmetry C(m,n)=C(n,m)\n";
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

    const long double answer_log10 = compute_log10_spanning_trees(options.m, options.n);
    std::cout << format_scientific_from_log10(answer_log10, 5) << '\n';
    return 0;
}
