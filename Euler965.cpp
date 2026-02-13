#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>

namespace {

long double compute_F(int N) {
    std::int64_t a = 0;
    std::int64_t b = 1;
    std::int64_t c = 1;
    std::int64_t d = N;

    long double sum = 0.0L;
    long double compensation = 0.0L;

    // On each Farey interval (a/b, c/d), f_N(x) = b*x - a.
    while (c <= N) {
        const long double term = 1.0L / (2.0L * static_cast<long double>(b) *
                                         static_cast<long double>(d) * static_cast<long double>(d));

        const long double y = term - compensation;
        const long double t = sum + y;
        compensation = (t - sum) - y;
        sum = t;

        const std::int64_t k = (N + b) / d;
        const std::int64_t next_c = k * c - a;
        const std::int64_t next_d = k * d - b;
        a = c;
        b = d;
        c = next_c;
        d = next_d;
    }

    return sum;
}

bool nearly_equal(long double x, long double y, long double eps) {
    return std::fabsl(x - y) <= eps;
}

bool run_checkpoints() {
    if (!nearly_equal(compute_F(1), 0.5L, 1e-18L)) {
        std::cerr << "Checkpoint failed for F(1)" << '\n';
        return false;
    }

    if (!nearly_equal(compute_F(4), 0.25L, 1e-18L)) {
        std::cerr << "Checkpoint failed for F(4)" << '\n';
        return false;
    }

    const long double f10 = compute_F(10);
    if (!nearly_equal(f10, 19.0L / 144.0L, 1e-18L)) {
        std::cerr << "Checkpoint failed for F(10): got " << std::setprecision(20) << f10 << '\n';
        return false;
    }

    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    const long double answer = compute_F(10'000);
    std::cout << std::fixed << std::setprecision(13) << answer << '\n';
    return 0;
}
