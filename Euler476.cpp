#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct Options {
    int n = 1803;
    bool run_checkpoints = true;
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
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.n < 2) {
        std::cerr << "--n must be at least 2.\n";
        return false;
    }
    return true;
}

long double solve_average(const int n) {
    constexpr long double kPi = 3.141592653589793238462643383279502884L;
    long double total = 0.0L;
    std::uint64_t count = 0;

    for (int a = 1; a <= n; ++a) {
        for (int b = a; a + b <= n; ++b) {
            const int c_max = a + b - 1;
            for (int c = b; c <= c_max; ++c) {
                const long double p = static_cast<long double>(a + b + c);
                const long double x = p - 2.0L * static_cast<long double>(a);  // b+c-a
                const long double y = p - 2.0L * static_cast<long double>(b);  // a+c-b
                const long double z = p - 2.0L * static_cast<long double>(c);  // a+b-c

                // Inradius squared: r^2 = ((b+c-a)(a+c-b)(a+b-c)) / (4(a+b+c))
                const long double r2 = (x * y * z) / (4.0L * p);

                // For angle A (opposite side a): sin(A/2)^2 = ((s-b)(s-c))/(bc)
                const long double sin_half_a =
                    std::sqrt((y * z) / (4.0L * static_cast<long double>(b) * c));
                // For angle B (opposite side b): sin(B/2)^2 = ((s-a)(s-c))/(ac)
                const long double sin_half_b =
                    std::sqrt((x * z) / (4.0L * static_cast<long double>(a) * c));

                const long double q_a = (1.0L - sin_half_a) / (1.0L + sin_half_a);
                const long double q_b = (1.0L - sin_half_b) / (1.0L + sin_half_b);

                // After placing the incircle, each corner yields a geometric chain of circles.
                // The largest extra circle is always in corner A (smallest angle), and the third
                // circle is either in corner B or the next one in corner A.
                const long double extra2 = q_a * q_a;
                const long double extra3 = std::max(q_b * q_b, extra2 * extra2);
                const long double area = kPi * r2 * (1.0L + extra2 + extra3);

                total += area;
                ++count;
            }
        }
    }

    return total / static_cast<long double>(count);
}

bool nearly_equal(const long double a, const long double b, const long double eps) {
    return std::fabsl(a - b) <= eps;
}

bool run_checkpoints() {
    if (!nearly_equal(solve_average(2), 0.31998L, 5e-5L)) {
        std::cerr << "Checkpoint failed: S(2)\n";
        return false;
    }
    if (!nearly_equal(solve_average(5), 1.25899L, 5e-5L)) {
        std::cerr << "Checkpoint failed: S(5)\n";
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

    const long double ans = solve_average(options.n);
    std::cout << std::fixed << std::setprecision(5) << static_cast<double>(ans) << '\n';
    return 0;
}
