#include <array>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Real = boost::multiprecision::cpp_dec_float_100;

constexpr int kStateCount = 20;  // 10 digits * run length {1,2}

int state_index(const int digit, const int run_len) {
    return digit * 2 + (run_len - 1);
}

bool has_three_equal_consecutive(std::uint64_t n) {
    int prev1 = -1;
    int prev2 = -2;
    while (n > 0) {
        const int d = static_cast<int>(n % 10ULL);
        n /= 10ULL;
        if (d == prev1 && d == prev2) {
            return true;
        }
        prev2 = prev1;
        prev1 = d;
    }
    return false;
}

struct Automaton {
    std::array<std::vector<std::pair<int, int>>, kStateCount> trans;
};

Automaton build_automaton() {
    Automaton a;
    for (int d = 0; d <= 9; ++d) {
        for (int r = 1; r <= 2; ++r) {
            const int s = state_index(d, r);
            a.trans[static_cast<std::size_t>(s)].clear();
            for (int x = 0; x <= 9; ++x) {
                if (x == d) {
                    if (r == 2) {
                        continue;
                    }
                    a.trans[static_cast<std::size_t>(s)].push_back({x, state_index(d, 2)});
                } else {
                    a.trans[static_cast<std::size_t>(s)].push_back({x, state_index(x, 1)});
                }
            }
        }
    }
    return a;
}

std::vector<std::vector<Real>> solve_moments(const Automaton& automaton, const int max_m) {
    std::vector<std::vector<Real>> h(static_cast<std::size_t>(kStateCount),
                                     std::vector<Real>(static_cast<std::size_t>(max_m + 1), Real(0)));

    // Binomial coefficients.
    std::vector<std::vector<std::uint64_t>> binom(static_cast<std::size_t>(max_m + 1),
                                                  std::vector<std::uint64_t>(static_cast<std::size_t>(max_m + 1), 0));
    for (int n = 0; n <= max_m; ++n) {
        binom[static_cast<std::size_t>(n)][0] = 1;
        binom[static_cast<std::size_t>(n)][static_cast<std::size_t>(n)] = 1;
        for (int k = 1; k < n; ++k) {
            binom[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)] =
                binom[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(k - 1)] +
                binom[static_cast<std::size_t>(n - 1)][static_cast<std::size_t>(k)];
        }
    }

    std::array<std::array<Real, 32>, 10> digit_pow{};
    for (int d = 0; d <= 9; ++d) {
        digit_pow[static_cast<std::size_t>(d)][0] = Real(1);
        for (int p = 1; p < 32; ++p) {
            digit_pow[static_cast<std::size_t>(d)][static_cast<std::size_t>(p)] =
                digit_pow[static_cast<std::size_t>(d)][static_cast<std::size_t>(p - 1)] * Real(d);
        }
    }

    for (int m = 0; m <= max_m; ++m) {
        const Real c = pow(Real(10), -(m + 1));

        std::vector<std::vector<Real>> a(static_cast<std::size_t>(kStateCount),
                                         std::vector<Real>(static_cast<std::size_t>(kStateCount + 1), Real(0)));

        for (int s = 0; s < kStateCount; ++s) {
            a[static_cast<std::size_t>(s)][static_cast<std::size_t>(s)] = Real(1);

            for (const auto& [digit, to] : automaton.trans[static_cast<std::size_t>(s)]) {
                (void)digit;
                a[static_cast<std::size_t>(s)][static_cast<std::size_t>(to)] -= c;
            }

            Real rhs = (m == 0) ? Real(1) : Real(0);
            for (const auto& [digit, to] : automaton.trans[static_cast<std::size_t>(s)]) {
                for (int u = 0; u < m; ++u) {
                    rhs += c * Real(binom[static_cast<std::size_t>(m)][static_cast<std::size_t>(u)]) *
                           digit_pow[static_cast<std::size_t>(digit)][static_cast<std::size_t>(m - u)] *
                           h[static_cast<std::size_t>(to)][static_cast<std::size_t>(u)];
                }
            }
            a[static_cast<std::size_t>(s)][static_cast<std::size_t>(kStateCount)] = rhs;
        }

        // Gaussian elimination with partial pivoting.
        for (int col = 0; col < kStateCount; ++col) {
            int pivot = col;
            for (int r = col + 1; r < kStateCount; ++r) {
                if (abs(a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]) >
                    abs(a[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(col)])) {
                    pivot = r;
                }
            }
            std::swap(a[static_cast<std::size_t>(col)], a[static_cast<std::size_t>(pivot)]);

            const Real pv = a[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
            for (int j = col; j <= kStateCount; ++j) {
                a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)] /= pv;
            }

            for (int r = 0; r < kStateCount; ++r) {
                if (r == col) {
                    continue;
                }
                const Real f = a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)];
                if (abs(f) < Real("1e-70")) {
                    continue;
                }
                for (int j = col; j <= kStateCount; ++j) {
                    a[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -=
                        f * a[static_cast<std::size_t>(col)][static_cast<std::size_t>(j)];
                }
            }
        }

        for (int s = 0; s < kStateCount; ++s) {
            h[static_cast<std::size_t>(s)][static_cast<std::size_t>(m)] =
                a[static_cast<std::size_t>(s)][static_cast<std::size_t>(kStateCount)];
        }
    }

    return h;
}

Real estimate_series_value(const std::vector<std::vector<Real>>& h, const int max_m, const int prefix_digits) {
    const std::uint64_t low = 1;
    std::uint64_t high_small = 1;
    for (int i = 0; i < prefix_digits - 1; ++i) {
        high_small *= 10ULL;
    }
    const std::uint64_t high_prefix = high_small * 10ULL;

    Real total = 0;

    for (std::uint64_t n = low; n < high_small; ++n) {
        if (!has_three_equal_consecutive(n)) {
            total += Real(1) / Real(n);
        }
    }

    for (std::uint64_t p = high_small; p < high_prefix; ++p) {
        if (has_three_equal_consecutive(p)) {
            continue;
        }

        const int last = static_cast<int>(p % 10ULL);
        const int prev = static_cast<int>((p / 10ULL) % 10ULL);
        const int run = (last == prev) ? 2 : 1;
        const int s = state_index(last, run);

        const Real inv_p = Real(1) / Real(p);
        Real inv_power = inv_p;
        Real contrib = 0;
        int sign = 1;
        for (int m = 0; m <= max_m; ++m) {
            const Real term = h[static_cast<std::size_t>(s)][static_cast<std::size_t>(m)] * inv_power;
            contrib += (sign > 0) ? term : -term;
            inv_power *= inv_p;
            sign = -sign;
        }
        total += contrib;
    }

    return total;
}

bool run_checkpoints() {
    int omitted = 0;
    for (int n = 1; n <= 1200; ++n) {
        if (has_three_equal_consecutive(static_cast<std::uint64_t>(n))) {
            ++omitted;
        }
    }
    if (omitted != 20) {
        std::cerr << "Checkpoint failed: omitted terms up to 1200\n";
        return false;
    }

    const Automaton automaton = build_automaton();
    const int max_m = 14;
    const auto h = solve_moments(automaton, max_m);

    const Real s2 = estimate_series_value(h, max_m, 2);
    const Real s3 = estimate_series_value(h, max_m, 3);
    if (abs(s2 - s3) > Real("1e-15")) {
        std::cerr << "Checkpoint failed: prefix split mismatch\n";
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

    const Automaton automaton = build_automaton();
    const int max_m = 16;
    const auto h = solve_moments(automaton, max_m);
    const Real value = estimate_series_value(h, max_m, 3);

    std::cout << std::fixed << std::setprecision(10) << value << '\n';
    return 0;
}
