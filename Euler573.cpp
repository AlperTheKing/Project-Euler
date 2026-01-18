#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {

long double compute_expected(long long n) {
    if (n <= 1) {
        return 1.0L;
    }

    long double nn = static_cast<long double>(n);
    long double c = std::powl((nn - 1.0L) / nn, nn - 1.0L);  // k=1 contribution.
    long double total = 1.0L + c;  // k=n contributes exactly 1.

    for (long long k = 1; k <= n - 2; ++k) {
        long long m = n - k;
        // Recurrence: c_{k+1} = c_k * ((k+1)/k)^k * ((m-1)/m)^(m-1).
        long double log_r = static_cast<long double>(k) * std::log1pl(1.0L / k) +
                            static_cast<long double>(m - 1) * std::log1pl(-1.0L / m);
        c *= std::expl(log_r);
        total += c;
    }
    return total;
}

bool validate() {
    struct Check {
        long long n;
        long double expected;
    };

    const std::vector<Check> checks = {
        {3, 17.0L / 9.0L},
        {4, 2.21875L},
        {5, 2.5104L},
        {10, 3.66021568L},
    };

    const long double tol = 1e-12L;
    for (const auto& chk : checks) {
        long double got = compute_expected(chk.n);
        if (std::fabsl(got - chk.expected) > tol) {
            std::cerr << "Validation failed for n=" << chk.n << ": got "
                      << std::setprecision(17) << static_cast<double>(got)
                      << ", expected " << static_cast<double>(chk.expected) << "\n";
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!validate()) {
        return 1;
    }

    long long n = 1'000'000;
#if defined(_WIN32)
    bool stdin_is_tty = _isatty(_fileno(stdin));
#else
    bool stdin_is_tty = isatty(fileno(stdin));
#endif
    if (argc > 1) {
        n = std::stoll(argv[1]);
    } else if (!stdin_is_tty) {
        long long tmp = 0;
        if (std::cin >> tmp) {
            n = tmp;
        }
    }

    long double ans = compute_expected(n);
    std::cout << std::fixed << std::setprecision(4) << static_cast<double>(ans) << "\n";
    return 0;
}
