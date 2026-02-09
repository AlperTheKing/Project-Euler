#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

constexpr u64 kLimit = 1000000000000000ULL;

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

std::vector<int> compute_totients(const int n_max) {
    std::vector<int> phi(static_cast<std::size_t>(n_max + 1));
    for (int i = 0; i <= n_max; ++i) {
        phi[static_cast<std::size_t>(i)] = i;
    }

    for (int p = 2; p <= n_max; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (int multiple = p; multiple <= n_max; multiple += p) {
            phi[static_cast<std::size_t>(multiple)] -= phi[static_cast<std::size_t>(multiple)] / p;
        }
    }

    return phi;
}

cpp_int central_binomial(const int n) {
    cpp_int c = 1;
    for (int i = 1; i <= n; ++i) {
        c *= (n + i);
        c /= i;
    }
    return c;
}

int compute_n_upper_bound(const u64 limit) {
    const cpp_int limit_big = limit;
    int n = 1;
    int best = 0;
    while (true) {
        const cpp_int bound = central_binomial(n) / (2 * n);
        if (bound > limit_big) {
            break;
        }
        best = n;
        ++n;
    }
    return best;
}

int compute_m_upper_bound(const u64 limit) {
    const cpp_int limit_big = limit;

    cpp_int factorial = 1;  // 1! = 1
    int best = 1;

    for (int m = 2;; ++m) {
        if (factorial > limit_big) {
            break;
        }
        best = m;
        factorial *= m;  // factorial now becomes m!
    }

    return best;
}

std::vector<cpp_int> factorials_up_to(const int n_max) {
    std::vector<cpp_int> fact(static_cast<std::size_t>(n_max + 1));
    fact[0] = 1;
    for (int i = 1; i <= n_max; ++i) {
        fact[static_cast<std::size_t>(i)] = fact[static_cast<std::size_t>(i - 1)] * i;
    }
    return fact;
}

std::vector<std::vector<cpp_int>> precompute_fact_powers(
    const std::vector<cpp_int>& fact,
    const int n_max,
    const int m_max
) {
    std::vector<std::vector<cpp_int>> powers(
        static_cast<std::size_t>(n_max + 1),
        std::vector<cpp_int>(static_cast<std::size_t>(m_max + 1), cpp_int(1))
    );

    for (int n = 0; n <= n_max; ++n) {
        powers[static_cast<std::size_t>(n)][0] = 1;
        for (int m = 1; m <= m_max; ++m) {
            powers[static_cast<std::size_t>(n)][static_cast<std::size_t>(m)] =
                powers[static_cast<std::size_t>(n)][static_cast<std::size_t>(m - 1)] *
                fact[static_cast<std::size_t>(n)];
        }
    }

    return powers;
}

cpp_int f_value(
    const int m,
    const int n,
    const std::vector<int>& phi,
    const std::vector<cpp_int>& fact,
    const std::vector<std::vector<cpp_int>>& fact_powers
) {
    const int total = m * n;
    cpp_int sum = 0;

    for (int l = 1; l <= n; ++l) {
        if ((n % l) != 0) {
            continue;
        }

        const int slice = n / l;
        const cpp_int numerator = fact[static_cast<std::size_t>(total / l)];
        const cpp_int denominator = fact_powers[static_cast<std::size_t>(slice)][static_cast<std::size_t>(m)];
        const cpp_int multinomial = numerator / denominator;

        sum += static_cast<cpp_int>(phi[static_cast<std::size_t>(l)]) * multinomial;
    }

    return sum / total;
}

std::vector<int> rotate_left(const std::vector<int>& v, const int shift) {
    const int n = static_cast<int>(v.size());
    std::vector<int> out(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        out[static_cast<std::size_t>(i)] = v[static_cast<std::size_t>((i + shift) % n)];
    }
    return out;
}

u64 brute_count_necklaces(const int m, const int n) {
    std::vector<int> arrangement;
    arrangement.reserve(static_cast<std::size_t>(m * n));
    for (int color = 0; color < m; ++color) {
        for (int cnt = 0; cnt < n; ++cnt) {
            arrangement.push_back(color);
        }
    }

    std::set<std::vector<int>> canonical;

    std::sort(arrangement.begin(), arrangement.end());
    do {
        std::vector<int> best = arrangement;
        for (int shift = 1; shift < m * n; ++shift) {
            std::vector<int> r = rotate_left(arrangement, shift);
            if (r < best) {
                best = std::move(r);
            }
        }
        canonical.insert(std::move(best));
    } while (std::next_permutation(arrangement.begin(), arrangement.end()));

    return static_cast<u64>(canonical.size());
}

bool run_checkpoints(
    const std::vector<int>& phi,
    const std::vector<cpp_int>& fact,
    const std::vector<std::vector<cpp_int>>& fact_powers
) {
    struct Case {
        int m;
        int n;
        u64 expected;
    };

    const std::vector<Case> explicit_cases = {
        {2, 1, 1},
        {2, 2, 2},
        {3, 1, 2},
        {3, 2, 16},
    };

    for (const Case& c : explicit_cases) {
        const cpp_int value = f_value(c.m, c.n, phi, fact, fact_powers);
        if (value != c.expected) {
            std::cerr << "Explicit checkpoint failed for f(" << c.m << ',' << c.n
                      << "): got " << value << ", expected " << c.expected << '\n';
            return false;
        }
    }

    const std::vector<std::pair<int, int>> brute_cases = {
        {2, 3},
        {2, 4},
        {3, 2},
        {4, 2},
    };

    for (const auto& [m, n] : brute_cases) {
        const cpp_int formula = f_value(m, n, phi, fact, fact_powers);
        const u64 brute = brute_count_necklaces(m, n);
        if (formula != brute) {
            std::cerr << "Brute checkpoint failed for f(" << m << ',' << n
                      << "): formula=" << formula << ", brute=" << brute << '\n';
            return false;
        }
    }

    return true;
}

u64 solve() {
    const int n_upper = compute_n_upper_bound(kLimit);
    const int m_upper = compute_m_upper_bound(kLimit);

    const int max_factorial_needed = m_upper * n_upper;
    const std::vector<cpp_int> fact = factorials_up_to(max_factorial_needed);
    const std::vector<int> phi = compute_totients(n_upper);
    const std::vector<std::vector<cpp_int>> fact_powers = precompute_fact_powers(fact, n_upper, m_upper);

    if (!run_checkpoints(phi, fact, fact_powers)) {
        return 0;
    }

    const cpp_int limit_big = kLimit;
    u64 answer = 0;

    for (int m = 2; m <= m_upper; ++m) {
        for (int n = 1; n <= n_upper; ++n) {
            const cpp_int value = f_value(m, n, phi, fact, fact_powers);
            if (value <= limit_big) {
                answer += value.convert_to<u64>();
            }
        }
    }

    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const int n_upper = compute_n_upper_bound(kLimit);
    const int m_upper = compute_m_upper_bound(kLimit);
    const int max_factorial_needed = m_upper * n_upper;

    const std::vector<cpp_int> fact = factorials_up_to(max_factorial_needed);
    const std::vector<int> phi = compute_totients(n_upper);
    const std::vector<std::vector<cpp_int>> fact_powers = precompute_fact_powers(fact, n_upper, m_upper);

    if (options.run_checkpoints && !run_checkpoints(phi, fact, fact_powers)) {
        return 2;
    }

    const cpp_int limit_big = kLimit;
    u64 answer = 0;

    for (int m = 2; m <= m_upper; ++m) {
        for (int n = 1; n <= n_upper; ++n) {
            const cpp_int value = f_value(m, n, phi, fact, fact_powers);
            if (value <= limit_big) {
                answer += value.convert_to<u64>();
            }
        }
    }

    std::cout << answer << '\n';
    return 0;
}
