#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    int limit_n = 50000000;
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
        if (parse_int_after_prefix(arg, "--limit-n=", options.limit_n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit_n >= 2;
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
    for (int p = 2; static_cast<u64>(p) * p <= static_cast<u64>(limit); ++p) {
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

int integer_sqrt_int(const int x) {
    int r = static_cast<int>(std::sqrt(static_cast<double>(x)));
    while (static_cast<long long>(r + 1) * (r + 1) <= x) {
        ++r;
    }
    while (static_cast<long long>(r) * r > x) {
        --r;
    }
    return r;
}

int solve(const int limit_n) {
    const u64 max_t = 2ULL * static_cast<u64>(limit_n) * static_cast<u64>(limit_n) - 1ULL;
    const int max_prime = static_cast<int>(std::sqrt(static_cast<long double>(max_t)));
    const std::vector<int> primes = sieve_primes(max_prime);

    std::vector<std::uint8_t> composite(static_cast<std::size_t>(limit_n + 1), 0U);

    for (int p : primes) {
        if (p == 2) {
            continue;
        }
        const int mod8 = p & 7;
        if (mod8 != 1 && mod8 != 7) {
            continue;
        }

        const u64 inv2 = static_cast<u64>((p + 1) / 2);
        u64 root = 0;
        if (mod8 == 7) {
            const u64 sqrt2 = mod_pow(2ULL, static_cast<u64>((p + 1) / 4), static_cast<u64>(p));
            root = static_cast<u64>((static_cast<u128>(sqrt2) * inv2) % static_cast<u64>(p));
        } else {
            root = tonelli_shanks(inv2, static_cast<u64>(p));
            if (root == 0) {
                continue;
            }
        }
        const int r1 = static_cast<int>(root);
        const int r2 = (r1 == 0) ? 0 : (p - r1);

        int skip_n = -1;
        const int half = (p + 1) / 2;
        const int sr = integer_sqrt_int(half);
        if (sr * sr == half) {
            skip_n = sr;
        }

        const auto mark_progression = [&](const int residue) {
            int n = residue;
            if (n < 2) {
                const int delta = 2 - n;
                n += ((delta + p - 1) / p) * p;
            }
            for (; n <= limit_n; n += p) {
                if (n == skip_n) {
                    continue;
                }
                composite[static_cast<std::size_t>(n)] = 1U;
            }
        };

        mark_progression(r1);
        if (r2 != r1) {
            mark_progression(r2);
        }
    }

    int count = 0;
    for (int n = 2; n <= limit_n; ++n) {
        if (composite[static_cast<std::size_t>(n)] == 0U) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(8) != 6) {
        std::cerr << "Checkpoint failed for n <= 8 prefix sample" << '\n';
        return false;
    }
    if (solve(10000) != 2202) {
        std::cerr << "Checkpoint failed for n <= 10000 sample" << '\n';
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
    std::cout << solve(options.limit_n) << '\n';
    return 0;
}
