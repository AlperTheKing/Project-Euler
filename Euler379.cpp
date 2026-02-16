#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr u64 kDefaultN = 1'000'000'000'000ULL;
constexpr u64 kCheckpointN = 1'000'000ULL;
constexpr u64 kCheckpointExpected = 37'429'395ULL;

struct Options {
    u64 n = kDefaultN;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            options.n = parsed_u64;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) > 0ULL && (r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 icbrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while ((r + 1ULL) > 0ULL &&
           static_cast<unsigned __int128>(r + 1ULL) * static_cast<unsigned __int128>(r + 1ULL) *
                   static_cast<unsigned __int128>(r + 1ULL) <=
               static_cast<unsigned __int128>(x)) {
        ++r;
    }
    while (static_cast<unsigned __int128>(r) * static_cast<unsigned __int128>(r) *
               static_cast<unsigned __int128>(r) >
           static_cast<unsigned __int128>(x)) {
        --r;
    }
    return r;
}

i64 sq_sum_floor(const i64 n, const i64 x1, const i64 x2) {
    if (x1 > x2) {
        return 0;
    }

    i64 x = x2;
    i64 s = 0;
    i64 beta = n / (x + 1);
    i64 epsilon = n % (x + 1);
    i64 delta = (n / x) - beta;
    i64 gamma = beta - x * delta;

    while (x >= x1) {
        epsilon += gamma;
        if (epsilon >= x) {
            delta += 1;
            gamma -= x;
            epsilon -= x;
            if (epsilon >= x) {
                delta += 1;
                gamma -= x;
                epsilon -= x;
                if (epsilon >= x) {
                    break;
                }
            }
        } else if (epsilon < 0) {
            delta -= 1;
            gamma += x;
            epsilon += x;
        }

        gamma += 2 * delta;
        beta += delta;
        s += beta;
        --x;
    }

    if (x < x1) {
        return s;
    }

    epsilon = n % (x + 1);
    delta = (n / x) - beta;
    gamma = beta - x * delta;

    while (x >= x1) {
        epsilon += gamma;
        const i64 delta2 = epsilon / x;
        delta += delta2;
        epsilon -= x * delta2;
        gamma += 2 * delta - x * delta2;
        beta += delta;
        s += beta;
        --x;
    }

    while (x >= x1) {
        s += n / x;
        --x;
    }

    return s;
}

i64 summatory_d3(const i64 n) {
    const i64 cbrt_n = static_cast<i64>(icbrt_u64(static_cast<u64>(n)));
    i64 acc = 0;

    for (i64 z = 1; z <= cbrt_n; ++z) {
        const i64 nz = n / z;
        const i64 sqrt_nz = static_cast<i64>(isqrt_u64(static_cast<u64>(nz)));
        acc += 2 * sq_sum_floor(nz, z + 1, sqrt_nz) - sqrt_nz * sqrt_nz + nz / z;
    }

    acc *= 3;
    acc += cbrt_n * cbrt_n * cbrt_n;
    return acc;
}

std::vector<int> mobius_sieve(const int n) {
    std::vector<int> mu(static_cast<std::size_t>(n) + 1ULL, 0);
    std::vector<int> lp(static_cast<std::size_t>(n) + 1ULL, 0);
    std::vector<int> primes;

    mu[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (lp[static_cast<std::size_t>(i)] == 0) {
            lp[static_cast<std::size_t>(i)] = i;
            primes.push_back(i);
            mu[static_cast<std::size_t>(i)] = -1;
        }
        for (const int p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > static_cast<u64>(n)) {
                break;
            }
            lp[static_cast<std::size_t>(v)] = p;
            if (i % p == 0) {
                mu[static_cast<std::size_t>(v)] = 0;
                break;
            }
            mu[static_cast<std::size_t>(v)] = -mu[static_cast<std::size_t>(i)];
        }
    }

    return mu;
}

u64 gcd_u64(u64 a, u64 b) {
    while (b != 0ULL) {
        const u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u64 brute_f_from_lcm(const u64 n) {
    u64 count = 0ULL;
    for (u64 x = 1ULL; x <= n; ++x) {
        for (u64 y = x; y <= n; ++y) {
            const u64 g = gcd_u64(x, y);
            const u64 l = (x / g) * y;
            if (l == n) {
                ++count;
            }
        }
    }
    return count;
}

u64 formula_f_from_factorization(u64 n) {
    u64 divisor_count_n2 = 1ULL;
    for (u64 p = 2ULL; p * p <= n; ++p) {
        if (n % p != 0ULL) {
            continue;
        }

        int exponent = 0;
        while (n % p == 0ULL) {
            n /= p;
            ++exponent;
        }
        divisor_count_n2 *= static_cast<u64>(2 * exponent + 1);
    }
    if (n > 1ULL) {
        divisor_count_n2 *= 3ULL;
    }
    return (divisor_count_n2 + 1ULL) / 2ULL;
}

u64 summatory_g_via_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n) + 1ULL, 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        for (u64 j = static_cast<u64>(i); j <= static_cast<u64>(n); j += static_cast<u64>(i)) {
            if (spf[static_cast<std::size_t>(j)] == 0) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }

    u64 total = 0ULL;
    for (int x = 1; x <= n; ++x) {
        int value = x;
        u64 divisor_count_n2 = 1ULL;
        while (value > 1) {
            const int p = spf[static_cast<std::size_t>(value)];
            int exponent = 0;
            while (value % p == 0) {
                value /= p;
                ++exponent;
            }
            divisor_count_n2 *= static_cast<u64>(2 * exponent + 1);
        }
        total += (divisor_count_n2 + 1ULL) / 2ULL;
    }
    return total;
}

u64 solve(const u64 n) {
    const i64 I = static_cast<i64>(icbrt_u64(n));
    const i64 D = static_cast<i64>(isqrt_u64(n / static_cast<u64>(I)));

    i64 res = 0;
    std::vector<int> mu = mobius_sieve(static_cast<int>(D));
    std::vector<i64> mertens_small(static_cast<std::size_t>(D) + 1ULL, 0);

    for (i64 d = 1; d <= D; ++d) {
        const int mud = mu[static_cast<std::size_t>(d)];
        if (mud != 0) {
            const u64 div = static_cast<u64>(d) * static_cast<u64>(d);
            res += static_cast<i64>(mud) * summatory_d3(static_cast<i64>(n / div));
        }
        mertens_small[static_cast<std::size_t>(d)] =
            mertens_small[static_cast<std::size_t>(d - 1)] + static_cast<i64>(mud);
    }

    if (I >= 2) {
        std::vector<i64> small_diff(static_cast<std::size_t>(I), 0);
        for (i64 i = 1; i < I; ++i) {
            small_diff[static_cast<std::size_t>(i)] = summatory_d3(i);
        }

        res -= mertens_small[static_cast<std::size_t>(D)] *
               small_diff[static_cast<std::size_t>(I - 1)];

        for (i64 i = I - 1; i >= 2; --i) {
            small_diff[static_cast<std::size_t>(i)] -= small_diff[static_cast<std::size_t>(i - 1)];
        }

        std::vector<i64> mertens_big(static_cast<std::size_t>(I), 0);
        for (i64 i = I - 1; i >= 1; --i) {
            const u64 v = isqrt_u64(n / static_cast<u64>(i));
            const u64 v_sqrt = isqrt_u64(v);

            i64 m = 1 - static_cast<i64>(v) +
                    static_cast<i64>(v_sqrt) * mertens_small[static_cast<std::size_t>(v_sqrt)];

            for (u64 d = 2; d <= v_sqrt; ++d) {
                const u64 q = v / d;
                if (q <= static_cast<u64>(D)) {
                    m -= mertens_small[static_cast<std::size_t>(q)];
                } else {
                    const u64 idx = static_cast<u64>(i) * d * d;
                    m -= mertens_big[static_cast<std::size_t>(idx)];
                }
                m -= (mertens_small[static_cast<std::size_t>(d)] -
                      mertens_small[static_cast<std::size_t>(d - 1)]) *
                     static_cast<i64>(q);
            }

            mertens_big[static_cast<std::size_t>(i)] = m;
            res += small_diff[static_cast<std::size_t>(i)] * m;
        }
    }

    res += static_cast<i64>(n);
    return static_cast<u64>(res / 2);
}

bool run_checkpoints() {
    for (u64 n = 1ULL; n <= 200ULL; ++n) {
        const u64 brute = brute_f_from_lcm(n);
        const u64 formula = formula_f_from_factorization(n);
        if (brute != formula) {
            std::cerr << "Identity check failed at n=" << n << ": brute f(n)=" << brute
                      << ", formula f(n)=" << formula << '\n';
            return false;
        }
    }

    struct Checkpoint {
        int n;
        u64 expected;
    };
    const std::vector<Checkpoint> small = {
        {10, 29ULL},
        {100, 647ULL},
        {1'000, 11'751ULL},
        {10'000, 186'991ULL},
    };

    for (const Checkpoint& cp : small) {
        const u64 got = summatory_g_via_spf(cp.n);
        if (got != cp.expected) {
            std::cerr << "Small checkpoint failed for g(" << cp.n << "): expected " << cp.expected
                      << ", got " << got << '\n';
            return false;
        }
    }

    const u64 got = solve(kCheckpointN);
    if (got != kCheckpointExpected) {
        std::cerr << "Main checkpoint failed for g(" << kCheckpointN << "): expected "
                  << kCheckpointExpected << ", got " << got << '\n';
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
        return 1;
    }

    std::cout << solve(options.n) << '\n';
    return 0;
}
