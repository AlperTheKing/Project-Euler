#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 n = 1'000'000'000'000ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) return false;
    u64 parsed = 0ULL;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') return false;
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) continue;
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL;
}

u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / (x + 1ULL)) ++x;
    while (x > 0ULL && x > n / x) --x;
    return x;
}

u64 icbrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::cbrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / ((x + 1ULL) * (x + 1ULL))) ++x;
    while (x > 0ULL && x > n / (x * x)) --x;
    return x;
}

bool pow_leq(u64 base, int exp, u64 limit) {
    u128 v = 1;
    for (int i = 0; i < exp; ++i) {
        v *= base;
        if (v > static_cast<u128>(limit)) return false;
    }
    return true;
}

u64 iroot7_u64(u64 n) {
    long double x = static_cast<long double>(n);
    u64 r = static_cast<u64>(std::pow(x, 1.0L / 7.0L));
    while (pow_leq(r + 1ULL, 7, n)) ++r;
    while (r > 0ULL && !pow_leq(r, 7, n)) --r;
    return r;
}

std::vector<int> prime_sieve(u64 n) {
    if (n < 2ULL) return {};
    const std::size_t m = static_cast<std::size_t>(n + 1ULL);
    std::vector<unsigned char> is_prime(m, 1);
    is_prime[0] = 0;
    is_prime[1] = 0;
    const u64 r = isqrt_u64(n);
    for (u64 p = 2; p <= r; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) continue;
        for (u64 q = p * p; q <= n; q += p) is_prime[static_cast<std::size_t>(q)] = 0;
    }
    std::vector<int> primes;
    primes.reserve(m / 10);
    for (u64 i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(static_cast<int>(i));
    }
    return primes;
}

struct PiTables {
    u64 n = 0;
    u64 v = 0;
    std::vector<u64> smalls;
    std::vector<u64> larges;
};

PiTables tabulate_pis(u64 n, const std::vector<int>& primes) {
    PiTables t;
    t.n = n;
    t.v = isqrt_u64(n);
    t.smalls.assign(static_cast<std::size_t>(t.v + 1ULL), 0ULL);
    t.larges.assign(static_cast<std::size_t>(t.v + 1ULL), 0ULL);

    for (u64 i = 1; i <= t.v; ++i) {
        t.smalls[static_cast<std::size_t>(i)] = i - 1ULL;
        t.larges[static_cast<std::size_t>(i)] = n / i - 1ULL;
    }

    for (int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        if (p > t.v) break;
        const u64 p_cnt = t.smalls[static_cast<std::size_t>(p - 1ULL)];
        const u64 q = p * p;
        if (q > n) break;
        const u64 end = std::min<u64>(t.v, n / q);

        for (u64 i = 1; i <= end; ++i) {
            const u64 d = i * p;
            if (d <= t.v) {
                t.larges[static_cast<std::size_t>(i)] -=
                    t.larges[static_cast<std::size_t>(d)] - p_cnt;
            } else {
                t.larges[static_cast<std::size_t>(i)] -=
                    t.smalls[static_cast<std::size_t>(n / d)] - p_cnt;
            }
        }

        for (u64 i = t.v; i >= q; --i) {
            t.smalls[static_cast<std::size_t>(i)] -=
                t.smalls[static_cast<std::size_t>(i / p)] - p_cnt;
            if (i == q) break;
        }
    }

    return t;
}

u64 solve(u64 n) {
    const u64 v = isqrt_u64(n);
    const std::vector<int> primes = prime_sieve(v);
    const PiTables pi = tabulate_pis(n, primes);

    u64 ans = 0ULL;
    const u64 cbrt_n = icbrt_u64(n);
    const u64 pi_cbrt_n = pi.smalls[static_cast<std::size_t>(cbrt_n)];

    for (u64 pi_idx = 0; pi_idx < pi_cbrt_n; ++pi_idx) {
        const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(pi_idx)]);
        if (p * p * p >= n) break;
        const u64 m = n / p;
        const u64 q_lim = pi.smalls[static_cast<std::size_t>(isqrt_u64(m))];

        for (u64 pj_idx = pi_idx + 1ULL; pj_idx < q_lim; ++pj_idx) {
            const u64 q = static_cast<u64>(primes[static_cast<std::size_t>(pj_idx)]);
            const u64 r = m / q;
            const u64 pi_r = (r <= pi.v)
                                 ? pi.smalls[static_cast<std::size_t>(r)]
                                 : pi.larges[static_cast<std::size_t>(p * q)];
            ans += pi_r - pi.smalls[static_cast<std::size_t>(q)];
        }
    }

    for (u64 pi_idx = 0; pi_idx < pi_cbrt_n; ++pi_idx) {
        const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(pi_idx)]);
        const u64 p3 = p * p * p;
        const u64 r = n / p3;
        if (r <= 1ULL) break;
        ans += (r <= pi.v)
                   ? pi.smalls[static_cast<std::size_t>(r)]
                   : pi.larges[static_cast<std::size_t>(p3)];
    }

    const u64 root4 = isqrt_u64(isqrt_u64(n));
    ans -= pi.smalls[static_cast<std::size_t>(root4)];
    ans += pi.smalls[static_cast<std::size_t>(iroot7_u64(n))];

    return ans;
}

u64 brute_count(u64 n) {
    u64 count = 0ULL;
    for (u64 x = 1ULL; x <= n; ++x) {
        u64 y = x;
        int d = 1;
        for (u64 p = 2ULL; p * p <= y; ++p) {
            if (y % p != 0ULL) continue;
            int e = 0;
            while (y % p == 0ULL) {
                y /= p;
                ++e;
            }
            d *= (e + 1);
        }
        if (y > 1ULL) d *= 2;
        if (d == 8) ++count;
    }
    return count;
}

bool run_checkpoints() {
    if (solve(100ULL) != 10ULL) {
        std::cerr << "Checkpoint failed: f(100)=10" << '\n';
        return false;
    }
    if (solve(1000ULL) != 180ULL) {
        std::cerr << "Checkpoint failed: f(1000)=180" << '\n';
        return false;
    }
    if (solve(1'000'000ULL) != 224'427ULL) {
        std::cerr << "Checkpoint failed: f(1e6)=224427" << '\n';
        return false;
    }
    if (solve(20'000ULL) != brute_count(20'000ULL)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for n=20000" << '\n';
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) return 1;
    if (options.run_checkpoints && !run_checkpoints()) return 2;
    std::cout << solve(options.n) << '\n';
    return 0;
}
