#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;

constexpr u64 MOD = 1'000'000'007ULL;

struct Options {
    int m = 200;
    int threads = 0;  // 0 => auto
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
    int parsed = 0;
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(ch - '0');
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
        if (parse_int_after_prefix(arg, "--m=", options.m) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.m >= 1 && options.threads >= 0;
}

u64 gcd_u64(u64 a, u64 b) {
    while (b != 0ULL) {
        const u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return result;
}

i64 extended_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = extended_gcd(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return g;
}

u64 mod_inverse(u64 a, const u64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    if (g != 1) {
        return 0ULL;
    }
    i64 v = x % static_cast<i64>(mod);
    if (v < 0) {
        v += static_cast<i64>(mod);
    }
    return static_cast<u64>(v);
}

class PrimeFactorizer {
public:
    explicit PrimeFactorizer(const int max_n) {
        spf_.resize(static_cast<std::size_t>(max_n + 1));
        for (int i = 0; i <= max_n; ++i) {
            spf_[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(i);
        }
        if (max_n >= 1) {
            spf_[1] = 1U;
        }
        for (int i = 2; i <= max_n / i; ++i) {
            if (spf_[static_cast<std::size_t>(i)] != static_cast<std::uint32_t>(i)) {
                continue;
            }
            for (int j = i * i; j <= max_n; j += i) {
                if (spf_[static_cast<std::size_t>(j)] == static_cast<std::uint32_t>(j)) {
                    spf_[static_cast<std::size_t>(j)] = static_cast<std::uint32_t>(i);
                }
            }
        }
    }

    u64 euler_phi(u64 n) const {
        u64 result = n;
        u64 x = n;
        while (x > 1ULL) {
            const u64 p = static_cast<u64>(spf_[static_cast<std::size_t>(x)]);
            while (x % p == 0ULL) {
                x /= p;
            }
            result -= result / p;
        }
        return result;
    }

    std::vector<u64> prime_factors_with_multiplicity(u64 n) const {
        std::vector<u64> factors;
        while (n > 1ULL) {
            const u64 p = static_cast<u64>(spf_[static_cast<std::size_t>(n)]);
            factors.push_back(p);
            n /= p;
        }
        return factors;
    }

private:
    std::vector<std::uint32_t> spf_;
};

u64 multiplicative_order(const u64 a, const u64 m, const PrimeFactorizer& factorizer) {
    if (gcd_u64(a, m) != 1ULL) {
        return 0ULL;
    }
    if (m == 1ULL) {
        return 1ULL;
    }

    u64 order = factorizer.euler_phi(m);
    std::vector<u64> factors = factorizer.prime_factors_with_multiplicity(order);
    std::sort(factors.begin(), factors.end());
    factors.erase(std::unique(factors.begin(), factors.end()), factors.end());

    for (const u64 p : factors) {
        while (order % p == 0ULL && mod_pow(a, order / p, m) == 1ULL) {
            order /= p;
        }
    }
    return order;
}

u64 min_l_magnitude(const u64 p, const u64 q) {
    const u64 limit = (q + p - 1ULL) / p;
    u64 pow10 = 1ULL;
    u64 l = 1ULL;
    while (pow10 < limit) {
        pow10 *= 10ULL;
        ++l;
    }
    return l;
}

u64 solve_pair(const u64 u, const u64 v, const PrimeFactorizer& factorizer) {
    const u64 p = u * u * u;
    const u64 q = v * v * v;
    if (p >= 10ULL * q) {
        return 0ULL;
    }

    const u64 k_val = 10ULL * q - p;
    const u64 min_l = min_l_magnitude(p, q);
    u64 best_l = std::numeric_limits<u64>::max();
    u64 best_d = 10ULL;

    for (u64 d = 1ULL; d <= 9ULL; ++d) {
        if (p * (d + 1ULL) >= 10ULL * q) {
            continue;
        }
        const u64 g = gcd_u64(k_val, d);
        const u64 m_prime = k_val / g;
        if (gcd_u64(10ULL, m_prime) != 1ULL) {
            continue;
        }
        const u64 l_base = multiplicative_order(10ULL, m_prime, factorizer);
        if (l_base == 0ULL) {
            continue;
        }
        const u64 l_cand = ((min_l + l_base - 1ULL) / l_base) * l_base;
        if (l_cand < best_l || (l_cand == best_l && d < best_d)) {
            best_l = l_cand;
            best_d = d;
        }
    }

    if (best_l == std::numeric_limits<u64>::max()) {
        return 0ULL;
    }

    const u64 ten_pow_l = mod_pow(10ULL, best_l, MOD);
    const u64 term1 = static_cast<u64>((static_cast<u128>(best_d % MOD) * (q % MOD)) % MOD);
    const u64 term2 = (ten_pow_l + MOD - 1ULL) % MOD;
    const u64 k_inv = mod_inverse(k_val % MOD, MOD);
    u64 out = static_cast<u64>((static_cast<u128>(term1) * term2) % MOD);
    out = static_cast<u64>((static_cast<u128>(out) * k_inv) % MOD);
    return out;
}

u64 compute_total_serial(const int m, const PrimeFactorizer& factorizer) {
    u64 total = 0ULL;
    for (int u = 1; u <= m; ++u) {
        const u64 uu = static_cast<u64>(u);
        for (int v = 1; v <= m; ++v) {
            if (gcd_u64(uu, static_cast<u64>(v)) == 1ULL) {
                total += solve_pair(uu, static_cast<u64>(v), factorizer);
                if (total >= MOD) {
                    total %= MOD;
                }
            }
        }
    }
    return total % MOD;
}

u64 compute_total_parallel(const int m, const PrimeFactorizer& factorizer, int threads) {
    if (threads <= 1) {
        return compute_total_serial(m, factorizer);
    }
    threads = std::min(threads, std::max(1, m));
    const int chunk = (m + threads - 1) / threads;
    std::vector<std::thread> workers;
    std::vector<u64> partial(static_cast<std::size_t>(threads), 0ULL);

    for (int t = 0; t < threads; ++t) {
        const int start = t * chunk + 1;
        if (start > m) {
            break;
        }
        const int end = std::min(m, (t + 1) * chunk);
        workers.emplace_back([&, t, start, end]() {
            u64 local = 0ULL;
            for (int u = start; u <= end; ++u) {
                const u64 uu = static_cast<u64>(u);
                for (int v = 1; v <= m; ++v) {
                    if (gcd_u64(uu, static_cast<u64>(v)) == 1ULL) {
                        local += solve_pair(uu, static_cast<u64>(v), factorizer);
                        if (local >= MOD) {
                            local %= MOD;
                        }
                    }
                }
            }
            partial[static_cast<std::size_t>(t)] = local % MOD;
        });
    }
    for (std::thread& th : workers) {
        th.join();
    }
    u64 total = 0ULL;
    for (u64 v : partial) {
        total += v;
        if (total >= MOD) {
            total %= MOD;
        }
    }
    return total % MOD;
}

bool run_checkpoints(const PrimeFactorizer& factorizer, const int threads) {
    if (compute_total_serial(1, factorizer) != 1ULL) {
        std::cerr << "Checkpoint failed: T(1)=1" << '\n';
        return false;
    }
    if (compute_total_serial(3, factorizer) != 262'429'173ULL) {
        std::cerr << "Checkpoint failed: T(3)=262429173" << '\n';
        return false;
    }
    if (threads > 1) {
        const u64 serial = compute_total_serial(6, factorizer);
        const u64 parallel = compute_total_parallel(6, factorizer, threads);
        if (serial != parallel) {
            std::cerr << "Checkpoint failed: deterministic serial/parallel equality" << '\n';
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

    const int threads = (options.threads > 0)
                            ? options.threads
                            : static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
    const u64 max_n = 10ULL * static_cast<u64>(options.m) * static_cast<u64>(options.m) *
                      static_cast<u64>(options.m);
    if (max_n > static_cast<u64>(std::numeric_limits<int>::max())) {
        std::cerr << "Requested m is too large for SPF table bounds" << '\n';
        return 1;
    }

    PrimeFactorizer factorizer(static_cast<int>(max_n));
    if (options.run_checkpoints && !run_checkpoints(factorizer, threads)) {
        return 2;
    }

    std::cout << compute_total_parallel(options.m, factorizer, threads) << '\n';
    return 0;
}
