#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

long double T(const int s, const long double p) {
    if (s <= 0) {
        return 0.0L;
    }
    if (s == 1) {
        return 1.0L;
    }

    const long double q = 1.0L - p;

    std::vector<long double> qpow(static_cast<std::size_t>(s + 1), 1.0L);
    for (int i = 1; i <= s; ++i) {
        qpow[static_cast<std::size_t>(i)] = qpow[static_cast<std::size_t>(i - 1)] * q;
    }

    std::vector<long double> A(static_cast<std::size_t>(s + 1), 0.0L);
    std::vector<long double> B(static_cast<std::size_t>(s + 1), 0.0L);

    A[0] = 0.0L;
    B[0] = 0.0L;
    A[1] = 1.0L;
    B[1] = 0.0L;  // known infected for size 1 => no test needed

    for (int n = 2; n <= s; ++n) {
        const long double pe = 1.0L - qpow[static_cast<std::size_t>(n)];  // P(at least one infected in n)

        // State B(n): at least one infected is guaranteed.
        long double bestB = 1e100L;
        for (int k = 1; k < n; ++k) {
            const long double pa = qpow[static_cast<std::size_t>(k)] *
                                   (1.0L - qpow[static_cast<std::size_t>(n - k)]) / pe;
            const long double pb = 1.0L - pa;

            const long double cost = 1.0L +
                                     pa * B[static_cast<std::size_t>(n - k)] +
                                     pb * (B[static_cast<std::size_t>(k)] + A[static_cast<std::size_t>(n - k)]);
            if (cost < bestB) {
                bestB = cost;
            }
        }
        B[static_cast<std::size_t>(n)] = bestB;

        // State A(n): no prior condition.
        long double bestA = 1.0L + (1.0L - qpow[static_cast<std::size_t>(n)]) * B[static_cast<std::size_t>(n)];
        for (int k = 1; k < n; ++k) {
            const long double cost = 1.0L +
                                     A[static_cast<std::size_t>(n - k)] +
                                     (1.0L - qpow[static_cast<std::size_t>(k)]) * B[static_cast<std::size_t>(k)];
            if (cost < bestA) {
                bestA = cost;
            }
        }
        A[static_cast<std::size_t>(n)] = bestA;
    }

    return A[static_cast<std::size_t>(s)];
}

long double solve() {
    long double total = 0.0L;
    for (int i = 1; i <= 50; ++i) {
        const long double p = static_cast<long double>(i) / 100.0L;
        total += T(10000, p);
    }
    return total;
}

bool approx_equal_6dp(const long double x, const long double expected) {
    return std::fabsl(x - expected) < 0.0000005L;
}

bool run_checkpoints() {
    if (!approx_equal_6dp(T(25, 0.02L), 4.155452L)) {
        std::cerr << "Checkpoint failed: T(25,0.02)\n";
        return false;
    }
    if (!approx_equal_6dp(T(25, 0.10L), 12.702124L)) {
        std::cerr << "Checkpoint failed: T(25,0.10)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const long double answer = solve();
    std::cout << std::fixed << std::setprecision(6) << answer << '\n';
    return 0;
}
