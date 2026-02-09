#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;

struct Options {
    int n = 10000000;
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
    try {
        value = std::stoi(tail);
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
    return options.n >= 2;
}

int isqrt_i64(i64 x) {
    int r = static_cast<int>(std::sqrt(static_cast<long double>(x)));
    while (1LL * (r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (1LL * r * r > x) {
        --r;
    }
    return r;
}

u64 brute_count_small(const int n) {
    std::unordered_map<u64, int> multiplicity;
    const int max_diag = isqrt_i64(n - 1);

    for (int a = 1; a <= max_diag; ++a) {
        for (int d = 1; d <= max_diag; ++d) {
            const int rem = (n - (a * a + d * d) - 1) / 2;
            if (rem < 1) {
                continue;
            }
            const int t = a + d;
            for (int b = 1; b <= rem; ++b) {
                const int c_max = rem / b;
                for (int c = 1; c <= c_max; ++c) {
                    const int a11 = a * a + b * c;
                    const int a22 = d * d + b * c;
                    if (a11 + a22 >= n) {
                        continue;
                    }
                    const int a12 = b * t;
                    const int a21 = c * t;
                    const u64 key = (static_cast<u64>(a11) << 48) |
                                    (static_cast<u64>(a12) << 32) |
                                    (static_cast<u64>(a21) << 16) |
                                    static_cast<u64>(a22);
                    ++multiplicity[key];
                }
            }
        }
    }

    u64 count = 0;
    for (const auto& kv : multiplicity) {
        if (kv.second >= 2) {
            ++count;
        }
    }
    return count;
}

u64 F_value(const int n) {
    const i64 two_n = 2LL * n;
    const int p_max = isqrt_i64(two_n - 1);
    const int n_global = isqrt_i64((two_n - 1) / 5);
    const int l_max = (n_global * n_global) / 4;

    std::vector<int> tau(static_cast<std::size_t>(l_max + 1), 0);
    for (int d = 1; d <= l_max; ++d) {
        for (int x = d; x <= l_max; x += d) {
            ++tau[static_cast<std::size_t>(x)];
        }
    }

    u64 total = 0;

    for (int p = 2; p <= p_max; ++p) {
        const int p2 = p * p;
        const bool p_odd = (p & 1) != 0;

        for (int q = 1; q < p; ++q) {
            if (std::gcd(p, q) != 1) {
                continue;
            }

            const int pq2 = p2 + q * q;
            const int n_max = isqrt_i64((two_n - 1) / pq2);
            const bool mixed_parity = ((p ^ q) & 1) != 0;

            for (int nn = 1; nn <= n_max; ++nn) {
                if (mixed_parity && ((nn & 1) != 0)) {
                    continue;
                }

                const int m_max = (q * nn - 1) / p;  // |m| < qn/p ensures a,d > 0.
                if (m_max < 0) {
                    continue;
                }

                const int parity = (p_odd && ((q & 1) != 0)) ? (nn & 1) : 0;
                int m = -m_max;
                if ((m & 1) != parity) {
                    ++m;
                }

                const int n2 = nn * nn;
                for (; m <= m_max; m += 2) {
                    const int l_num = n2 - m * m;
                    if (l_num <= 0) {
                        continue;
                    }
                    const int l = l_num / 4;
                    total += static_cast<u64>(tau[static_cast<std::size_t>(l)]);
                }
            }
        }
    }

    return total;
}

bool run_checkpoints() {
    if (F_value(50) != 7ULL) {
        std::cerr << "Checkpoint failed: F(50)\n";
        return false;
    }
    if (F_value(1000) != 1019ULL) {
        std::cerr << "Checkpoint failed: F(1000)\n";
        return false;
    }
    if (F_value(200) != brute_count_small(200)) {
        std::cerr << "Checkpoint failed: formula vs brute at N=200\n";
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

    std::cout << F_value(options.n) << '\n';
    return 0;
}
