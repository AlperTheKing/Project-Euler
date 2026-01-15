#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

using int64 = long long;
using i128 = __int128_t;

struct Group {
    int64 l;
    int64 r;
    int64 k;
};

static int64 floor_sqrt(int64 x) {
    if (x <= 0) return 0;
    int64 r = static_cast<int64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) > 0 && (r + 1) * (r + 1) <= x) ++r;
    while (r * r > x) --r;
    return r;
}

static int64 count_gaussian_points(int64 x) {
    if (x < 0) return 0;
    const int64 r = floor_sqrt(x);
    int64 b = r;
    int64 b2 = b * b;
    int64 a2 = 0;
    int64 sum = 0;
    for (int64 a = 0; a <= r; ++a) {
        while (a2 + b2 > x) {
            --b;
            b2 -= 2 * b + 1;
        }
        sum += b;
        a2 += 2 * a + 1;
    }
    return 1 + 4 * sum;
}

static std::vector<int> build_spf(int limit) {
    std::vector<int> spf(limit + 1, 0);
    std::vector<int> primes;
    primes.reserve(limit / 10);
    for (int i = 2; i <= limit; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            if (p > spf[i] || i * 1LL * p > limit) break;
            spf[i * p] = p;
        }
    }
    return spf;
}

// Möbius contribution for proper Gaussian integers at a rational prime power.
static int64 local_g(int p, int exp) {
    if (p == 2) {
        if (exp == 1) return -1;
        return 0;
    }
    if (p % 4 == 1) {
        if (exp == 1) return -2;
        if (exp == 2) return 1;
        return 0;
    }
    if (exp == 2) return -1;
    return 0;
}

static int64 compute_f(int64 n, unsigned threads) {
    const int64 limit = floor_sqrt(n);
    std::vector<int> spf = build_spf(static_cast<int>(limit));

    std::vector<int64> g(limit + 1, 0);
    g[1] = 1;
    for (int64 m = 2; m <= limit; ++m) {
        int p = spf[m];
        int64 t = m;
        int exp = 0;
        while (t % p == 0) {
            t /= p;
            ++exp;
        }
        int64 local = local_g(p, exp);
        if (local == 0 || g[t] == 0) {
            g[m] = 0;
        } else {
            g[m] = local * g[t];
        }
    }

    std::vector<int64> prefix(limit + 1, 0);
    for (int64 i = 1; i <= limit; ++i) {
        prefix[i] = prefix[i - 1] + g[i];
    }

    std::vector<Group> groups;
    groups.reserve(static_cast<size_t>(2 * std::cbrt(static_cast<long double>(n)) + 10));
    for (int64 m = 1; m <= limit;) {
        int64 k = n / (m * m);
        int64 max_m = floor_sqrt(n / k);
        groups.push_back({m, max_m, k});
        m = max_m + 1;
    }

    if (threads == 0) threads = 1;
    if (threads > groups.size()) threads = static_cast<unsigned>(groups.size());
    if (threads == 0) threads = 1;

    std::vector<int64> A_values(groups.size(), 0);
    std::atomic<size_t> next_idx(0);
    auto worker = [&]() {
        while (true) {
            size_t idx = next_idx.fetch_add(1, std::memory_order_relaxed);
            if (idx >= groups.size()) break;
            A_values[idx] = count_gaussian_points(groups[idx].k);
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back(worker);
    }
    for (auto& th : pool) th.join();

    i128 total = 0;
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& gk = groups[i];
        int64 sum_g = prefix[gk.r] - prefix[gk.l - 1];
        i128 term = static_cast<i128>(sum_g) * (static_cast<i128>(A_values[i]) - 1);
        total += term;
    }

    // Divide by 4 to convert from all associates to proper representatives.
    return static_cast<int64>(total / 4);
}

static bool run_validations() {
    struct Check {
        int64 n;
        int64 expected;
    };
    const std::vector<Check> checks = {
        {10, 7},
        {100, 54},
        {10000, 5218},
        {100000000, 52126906},
    };

    for (const auto& check : checks) {
        int64 got = compute_f(check.n, 1);
        if (got != check.expected) {
            std::cerr << "Validation failed for f(" << check.n << "): got=" << got
                      << " expected=" << check.expected << "\n";
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int64 n = 100000000000000LL;
    if (argc > 1) n = std::atoll(argv[1]);
    if (n <= 0) {
        std::cerr << "n must be positive.\n";
        return 1;
    }

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 4;
    if (argc > 2) threads = static_cast<unsigned>(std::max(1, std::atoi(argv[2])));

    if (!run_validations()) return 1;

    int64 answer = compute_f(n, threads);
    std::cout << answer << "\n";
    return 0;
}
