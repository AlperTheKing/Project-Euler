#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
    u64 limit = 10'000'000'000ULL;
    bool run_checkpoints = true;
};

struct MEntry {
    int m = 0;
    u64 q = 0;
    u64 x_max = 0;
};

struct PellSolution {
    u64 x = 0;
    u64 y = 0;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    u64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(c - '0');
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
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1;
}

std::vector<int> build_primes(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    primes.reserve(n / 10);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(i) * p;
            if (v > n || p > spf[static_cast<std::size_t>(i)]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }
    return primes;
}

std::pair<u64, u64> squarefree_and_square_part(u64 n, const std::vector<int>& primes) {
    u64 s = 1;
    u64 q = 1;
    u64 x = n;

    for (int p : primes) {
        const u64 pp = static_cast<u64>(p);
        if (pp * pp > x) {
            break;
        }
        if (x % pp != 0) {
            continue;
        }
        int e = 0;
        while (x % pp == 0) {
            x /= pp;
            ++e;
        }
        if ((e & 1) != 0) {
            s *= pp;
        }
        for (int i = 0; i < e / 2; ++i) {
            q *= pp;
        }
    }

    if (x > 1) {
        s *= x;
    }

    return {s, q};
}

u64 floor_sqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (r * r > x) {
        --r;
    }
    return r;
}

bool find_fundamental_pell_limited(u64 s, u64 x_limit, u64& x1, u64& y1) {
    const u64 a0 = floor_sqrt_u64(s);
    if (a0 * a0 == s) {
        return false;
    }

    i128 m = 0;
    i128 d = 1;
    i128 a = static_cast<i128>(a0);

    i128 p_prev2 = 0;
    i128 p_prev1 = 1;
    i128 q_prev2 = 1;
    i128 q_prev1 = 0;

    for (int iter = 0; iter < 2'000'000; ++iter) {
        const i128 p = a * p_prev1 + p_prev2;
        const i128 q = a * q_prev1 + q_prev2;

        const i128 diff = p * p - static_cast<i128>(s) * q * q;
        if (diff == 1 && p > 1 && p <= static_cast<i128>(x_limit)) {
            x1 = static_cast<u64>(p);
            y1 = static_cast<u64>(q);
            return true;
        }

        if (p > static_cast<i128>(x_limit)) {
            return false;
        }

        p_prev2 = p_prev1;
        p_prev1 = p;
        q_prev2 = q_prev1;
        q_prev1 = q;

        m = d * a - m;
        d = (static_cast<i128>(s) - m * m) / d;
        a = (static_cast<i128>(a0) + m) / d;
    }

    return false;
}

std::vector<PellSolution> pell_solutions_up_to(u64 s, u64 x_limit) {
    u64 x1 = 0;
    u64 y1 = 0;
    if (!find_fundamental_pell_limited(s, x_limit, x1, y1)) {
        return {};
    }

    std::vector<PellSolution> out;
    u64 x = x1;
    u64 y = y1;

    while (x <= x_limit) {
        out.push_back({x, y});

        const i128 nx = static_cast<i128>(x1) * x + static_cast<i128>(s) * y1 * y;
        const i128 ny = static_cast<i128>(x1) * y + static_cast<i128>(y1) * x;
        if (nx > static_cast<i128>(x_limit)) {
            break;
        }
        x = static_cast<u64>(nx);
        y = static_cast<u64>(ny);
    }

    return out;
}

std::vector<u64> generate_square_pivots(const u64 limit) {
    const int m_max = static_cast<int>((std::sqrt(2.0L * static_cast<long double>(limit) + 1.0L) - 1.0L) / 2.0L);
    const std::vector<int> primes = build_primes(m_max + 2);

    std::unordered_map<u64, std::vector<MEntry>> grouped;
    grouped.reserve(static_cast<std::size_t>(m_max * 2));

    for (int m = 1; m <= m_max; ++m) {
        const u64 n = static_cast<u64>(m) * static_cast<u64>(m + 1);
        const auto [s, q] = squarefree_and_square_part(n, primes);
        const u64 x_max = (2ULL * limit - static_cast<u64>(m)) / static_cast<u64>(m);
        grouped[s].push_back({m, q, x_max});
    }

    std::vector<u64> pivots;
    pivots.reserve(1 << 20);

    for (auto& kv : grouped) {
        const u64 s = kv.first;
        std::vector<MEntry>& entries = kv.second;

        u64 x_max = 0;
        for (const MEntry& e : entries) {
            x_max = std::max(x_max, e.x_max);
        }

        const std::vector<PellSolution> sols = pell_solutions_up_to(s, x_max);
        if (sols.empty()) {
            continue;
        }

        for (const MEntry& e : entries) {
            const u64 min_d = static_cast<u64>(2 * e.m + 1);
            const i128 m128 = static_cast<i128>(e.m);
            const i128 qs128 = static_cast<i128>(e.q) * static_cast<i128>(s);

            for (const PellSolution& sol : sols) {
                if (sol.x < min_d) {
                    continue;
                }
                if ((sol.x & 1ULL) == 0ULL) {
                    continue;
                }

                const i128 numerator = m128 * static_cast<i128>(sol.x) +
                                       qs128 * static_cast<i128>(sol.y) + m128;
                if ((numerator & 1) != 0) {
                    continue;
                }

                const i128 k = numerator / 2;
                if (k > static_cast<i128>(limit)) {
                    continue;
                }
                if (k > 0) {
                    pivots.push_back(static_cast<u64>(k));
                }
            }
        }
    }

    std::sort(pivots.begin(), pivots.end());
    pivots.erase(std::unique(pivots.begin(), pivots.end()), pivots.end());
    return pivots;
}

u64 solve(const u64 limit) {
    const std::vector<u64> pivots = generate_square_pivots(limit);
    u64 sum = 0;
    for (u64 k : pivots) {
        sum += k;
    }
    return sum;
}

std::vector<u64> brute_pivots(const int limit) {
    std::vector<u64> pivots;
    for (int k = 1; k <= limit; ++k) {
        bool ok = false;
        for (int m = 1; m < k && !ok; ++m) {
            const i128 lhs = static_cast<i128>(m + 1) * k * (k - m);
            if (lhs % m != 0) {
                continue;
            }
            const i128 c = lhs / m;
            const i128 d = static_cast<i128>(m + 1) * (m + 1) + 4 * c;
            const i64 s = static_cast<i64>(std::sqrt(static_cast<long double>(d)));
            i64 root = s;
            while (static_cast<i128>(root + 1) * (root + 1) <= d) {
                ++root;
            }
            while (static_cast<i128>(root) * root > d) {
                --root;
            }
            if (static_cast<i128>(root) * root != d) {
                continue;
            }
            const i64 numer = -static_cast<i64>(m + 1) + root;
            if ((numer & 1) != 0) {
                continue;
            }
            const i64 n = numer / 2;
            if (n >= k) {
                ok = true;
            }
        }
        if (ok) {
            pivots.push_back(static_cast<u64>(k));
        }
    }
    return pivots;
}

bool run_checkpoints() {
    const auto fast_small = generate_square_pivots(5000);
    const auto brute_small = brute_pivots(5000);
    if (fast_small != brute_small) {
        std::cerr << "Checkpoint failed for limit=5000" << '\n';
        return false;
    }

    const auto sample = generate_square_pivots(120);
    const std::vector<u64> required = {4, 21, 24, 110};
    for (u64 x : required) {
        if (!std::binary_search(sample.begin(), sample.end(), x)) {
            std::cerr << "Checkpoint failed: missing sample pivot " << x << '\n';
            return false;
        }
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
    std::cout << solve(options.limit) << '\n';
    return 0;
}
