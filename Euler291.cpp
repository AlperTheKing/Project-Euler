#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 5000000000000000ULL;
    bool run_checkpoints = true;
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
    return options.limit >= 5ULL;
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        exp >>= 1U;
    }
    return result;
}

u64 tonelli_shanks(u64 n, u64 p) {
    if (n == 0) {
        return 0;
    }
    if (p == 2) {
        return n;
    }
    if (mod_pow(n, (p - 1) / 2, p) != 1) {
        return 0;
    }
    if (p % 4 == 3) {
        return mod_pow(n, (p + 1) / 4, p);
    }

    u64 q = p - 1;
    int s = 0;
    while ((q & 1ULL) == 0ULL) {
        q >>= 1U;
        ++s;
    }

    u64 z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) {
        ++z;
    }

    u64 m = static_cast<u64>(s);
    u64 c = mod_pow(z, q, p);
    u64 t = mod_pow(n, q, p);
    u64 r = mod_pow(n, (q + 1) / 2, p);

    while (t != 1) {
        u64 i = 1;
        u64 t2 = static_cast<u64>((static_cast<u128>(t) * t) % p);
        while (t2 != 1) {
            t2 = static_cast<u64>((static_cast<u128>(t2) * t2) % p);
            ++i;
        }
        const u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        r = static_cast<u64>((static_cast<u128>(r) * b) % p);
        t = static_cast<u64>((static_cast<u128>(t) * b % p * b) % p);
        c = static_cast<u64>((static_cast<u128>(b) * b) % p);
        m = i;
    }
    return r;
}

std::vector<int> sieve_primes(const int limit) {
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit + 1), 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;
    for (int p = 2; static_cast<u64>(p) * static_cast<u64>(p) <= static_cast<u64>(limit); ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (int q = p * p; q <= limit; q += p) {
            is_prime[static_cast<std::size_t>(q)] = 0U;
        }
    }
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::log(static_cast<double>(limit))));
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] != 0U) {
            primes.push_back(p);
        }
    }
    return primes;
}

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) > r && static_cast<u128>(r + 1) * static_cast<u128>(r + 1) <= x) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) > x) {
        --r;
    }
    return r;
}

u64 solve(const u64 limit) {
    // Prime form: p = 2n^2 + 2n + 1
    const u64 n_max = (isqrt_u64(2ULL * limit - 1ULL) - 1ULL) / 2ULL;
    const u64 max_q = 2ULL * n_max * n_max + 2ULL * n_max + 1ULL;
    const int max_prime = static_cast<int>(std::sqrt(static_cast<long double>(max_q)));
    const std::vector<int> primes = sieve_primes(max_prime);

    std::vector<std::uint8_t> composite(static_cast<std::size_t>(n_max + 1ULL), 0U);

    for (int p : primes) {
        if (p == 2) {
            continue;
        }
        if ((p & 3) != 1) {
            continue;  // -1 is non-residue.
        }

        const u64 sqrt_minus_one = tonelli_shanks(static_cast<u64>(p - 1), static_cast<u64>(p));
        if (sqrt_minus_one == 0ULL) {
            continue;
        }

        const u64 inv2 = static_cast<u64>((p + 1) / 2);
        const int r1 = static_cast<int>((static_cast<u128>(p - 1 + sqrt_minus_one) * inv2) % static_cast<u64>(p));
        const int r2 = static_cast<int>((static_cast<u128>(p - 1 + (p - sqrt_minus_one)) * inv2) %
                                        static_cast<u64>(p));

        const auto mark_progression = [&](const int residue) {
            u64 n = static_cast<u64>(residue);
            if (n == 0ULL) {
                n += static_cast<u64>(p);
            }
            for (; n <= n_max; n += static_cast<u64>(p)) {
                const u64 q = 2ULL * n * n + 2ULL * n + 1ULL;
                if (q == static_cast<u64>(p)) {
                    continue;  // prime itself, not composite.
                }
                composite[static_cast<std::size_t>(n)] = 1U;
            }
        };
        mark_progression(r1);
        if (r2 != r1) {
            mark_progression(r2);
        }
    }

    u64 count = 0;
    for (u64 n = 1; n <= n_max; ++n) {
        const u64 q = 2ULL * n * n + 2ULL * n + 1ULL;
        if (q >= limit) {
            break;
        }
        if (composite[static_cast<std::size_t>(n)] == 0U) {
            ++count;
        }
    }
    return count;
}

u64 brute_small(const u64 limit) {
    auto is_prime = [](u64 x) {
        if (x < 2ULL) {
            return false;
        }
        if ((x & 1ULL) == 0ULL) {
            return x == 2ULL;
        }
        for (u64 p = 3; p * p <= x; p += 2ULL) {
            if (x % p == 0ULL) {
                return false;
            }
        }
        return true;
    };

    u64 count = 0;
    for (u64 n = 1;; ++n) {
        const u64 q = 2ULL * n * n + 2ULL * n + 1ULL;
        if (q >= limit) {
            break;
        }
        if (is_prime(q)) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(1000ULL) != brute_small(1000ULL)) {
        std::cerr << "Checkpoint failed for limit=1000" << '\n';
        return false;
    }
    if (solve(1000000ULL) != brute_small(1000000ULL)) {
        std::cerr << "Checkpoint failed for limit=1,000,000" << '\n';
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
    std::cout << solve(options.limit) << '\n';
    return 0;
}
