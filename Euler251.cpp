#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using i128 = __int128_t;

struct Options {
    i64 limit = 110000000;
    int threads = 0;
    bool run_checkpoints = true;
};

bool parse_i64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    i64 parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<i64>(c - '0');
    }
    value = parsed;
    return true;
}

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    int parsed = 0;
    for (char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        parsed = parsed * 10 + static_cast<int>(c - '0');
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
        if (parse_i64_after_prefix(arg, "--limit=", options.limit) ||
            parse_int_after_prefix(arg, "--threads=", options.threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 1 && options.threads >= 0;
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    is_prime[0] = false;
    if (limit >= 1) {
        is_prime[1] = false;
    }

    for (int i = 2; static_cast<i64>(i) * i <= limit; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

std::vector<int> build_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    std::vector<int> primes;
    spf[0] = 0;
    if (n >= 1) {
        spf[1] = 1;
    }

    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(p) * static_cast<i64>(i);
            if (v > n || p > spf[static_cast<std::size_t>(i)]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }

    return spf;
}

void factorize_with_spf(int x, const std::vector<int>& spf, std::vector<std::pair<int, int>>& factors) {
    factors.clear();
    while (x > 1) {
        const int p = spf[static_cast<std::size_t>(x)];
        int e = 0;
        do {
            x /= p;
            ++e;
        } while (x > 1 && spf[static_cast<std::size_t>(x)] == p);
        factors.push_back({p, e});
    }
}

void factorize_square_part(i64 m, const std::vector<int>& primes, std::vector<std::pair<int, int>>& square_factors) {
    square_factors.clear();
    i64 x = m;

    for (int p : primes) {
        if (static_cast<i64>(p) * static_cast<i64>(p) > x) {
            break;
        }
        if ((x % p) != 0) {
            continue;
        }

        int e = 0;
        do {
            x /= p;
            ++e;
        } while ((x % p) == 0);

        if (e >= 2) {
            square_factors.push_back({p, e / 2});
        }
    }
}

void generate_divisors(const std::vector<std::pair<int, int>>& factors, std::vector<i64>& divisors) {
    divisors.clear();
    divisors.push_back(1);

    for (const auto& [p, e] : factors) {
        const std::size_t base_size = divisors.size();
        i64 pe = 1;
        for (int k = 1; k <= e; ++k) {
            pe *= static_cast<i64>(p);
            for (std::size_t i = 0; i < base_size; ++i) {
                divisors.push_back(divisors[i] * pe);
            }
        }
    }
}

i64 floor_sqrt_i64(const i64 x) {
    if (x <= 0) {
        return 0;
    }
    i64 r = static_cast<i64>(std::sqrt(static_cast<long double>(x)));
    while (static_cast<i128>(r + 1) * static_cast<i128>(r + 1) <= x) {
        ++r;
    }
    while (static_cast<i128>(r) * static_cast<i128>(r) > x) {
        --r;
    }
    return r;
}

bool can_have_solution(const i64 u, const i64 limit) {
    const long double uu = static_cast<long double>(u);
    const long double mm = 8.0L * uu - 3.0L;
    const long double lower_bound = 3.0L * uu - 1.0L + 3.0L * std::cbrt((uu * uu * mm) / 4.0L);
    return lower_bound <= static_cast<long double>(limit) + 1e-9L;
}

i64 max_u_bound(const i64 limit) {
    i64 lo = 0;
    i64 hi = (limit + 1) / 3;
    while (lo < hi) {
        const i64 mid = lo + (hi - lo + 1) / 2;
        if (can_have_solution(mid, limit)) {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    return lo;
}

i64 solve(const i64 limit, int threads) {
    if (limit < 8) {
        return 0;
    }

    const i64 max_u = max_u_bound(limit);
    if (max_u <= 0) {
        return 0;
    }

    const int sqrt_max_m = static_cast<int>(std::sqrt(static_cast<long double>(8 * max_u - 3))) + 1;
    const std::vector<int> primes = sieve_primes(sqrt_max_m);
    const std::vector<int> spf = build_spf(static_cast<int>(max_u));

    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
        if (threads <= 0) {
            threads = 1;
        }
    }
    threads = std::max(1, std::min(threads, static_cast<int>(max_u)));

    std::vector<i64> local_counts(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));

    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            std::vector<std::pair<int, int>> factors_u;
            std::vector<std::pair<int, int>> square_factors_m;
            std::vector<i64> divisors_h;
            std::vector<i64> divisors_b1;
            factors_u.reserve(16);
            square_factors_m.reserve(8);
            divisors_h.reserve(256);
            divisors_b1.reserve(64);

            i64 subtotal = 0;
            for (i64 u = static_cast<i64>(t) + 1; u <= max_u; u += threads) {
                const i64 a = 3 * u - 1;
                const i64 remaining = limit - a;
                if (remaining <= 1) {
                    continue;
                }

                factorize_with_spf(static_cast<int>(u), spf, factors_u);
                generate_divisors(factors_u, divisors_h);

                const i64 m = 8 * u - 3;
                factorize_square_part(m, primes, square_factors_m);
                generate_divisors(square_factors_m, divisors_b1);

                for (const i64 b1 : divisors_b1) {
                    const i64 b1_sq = b1 * b1;
                    const i64 m_over = m / b1_sq;

                    if (b1 >= remaining || b1 + m_over > remaining) {
                        continue;
                    }

                    const i64 h_max = floor_sqrt_i64((remaining - b1) / m_over);
                    if (h_max <= 0) {
                        continue;
                    }

                    i64 h_min = 1;
                    const i64 denom = remaining - m_over;
                    if (denom > 0) {
                        h_min = (u * b1 + denom - 1) / denom;
                    }

                    for (const i64 h : divisors_h) {
                        if (h < h_min || h > h_max) {
                            continue;
                        }
                        if (std::gcd(h, b1) != 1) {
                            continue;
                        }

                        const i64 b = (u / h) * b1;
                        const i64 c = h * h * m_over;
                        if (static_cast<i128>(b) + static_cast<i128>(c) <= remaining) {
                            ++subtotal;
                        }
                    }
                }
            }
            local_counts[static_cast<std::size_t>(t)] = subtotal;
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    i64 total = 0;
    for (const i64 v : local_counts) {
        total += v;
    }
    return total;
}

bool run_checkpoints() {
    if (solve(8, 1) != 1) {
        std::cerr << "Checkpoint failed for limit=8" << '\n';
        return false;
    }
    if (solve(100, 1) != 11) {
        std::cerr << "Checkpoint failed for limit=100" << '\n';
        return false;
    }
    if (solve(1000, 1) != 149) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
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

    std::cout << solve(options.limit, options.threads) << '\n';
    return 0;
}
