#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    int k = 10'000;
    u64 n = 1'000'000'000'000ULL;
    int p_lo = 2'000'000'000;
    int p_hi = 2'000'002'000;
    bool run_checkpoints = true;
};

bool parse_int_after_prefix(const std::string& arg, const std::string& prefix, int& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        out = static_cast<u64>(std::stoull(tail));
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
        if (parse_int_after_prefix(arg, "--k=", options.k)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--p-lo=", options.p_lo)) {
            continue;
        }
        if (parse_int_after_prefix(arg, "--p-hi=", options.p_hi)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    if (options.k < 0 || options.p_lo > options.p_hi || options.p_lo < 2) {
        std::cerr << "Invalid arguments.\n";
        return false;
    }
    return true;
}

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 out = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if (exp & 1ULL) {
            out = static_cast<u64>((static_cast<u128>(out) * base) % mod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % mod);
        exp >>= 1ULL;
    }
    return out;
}

u64 mod_inv(const u64 x, const u64 p) { return mod_pow(x, p - 2ULL, p); }

u64 add_mod(const u64 a, const u64 b, const u64 p) {
    const u64 c = a + b;
    return (c >= p || c < a) ? (c % p) : c;
}

u64 sub_mod(const u64 a, const u64 b, const u64 p) {
    return (a >= b) ? (a - b) : (a + p - b);
}

bool is_prime_64(const u64 n) {
    if (n < 2ULL) {
        return false;
    }
    for (const u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL,
                        31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0ULL) {
            return false;
        }
    }

    u64 d = n - 1ULL;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    auto witness = [&](const u64 a) {
        if (a % n == 0ULL) {
            return false;
        }
        u64 x = mod_pow(a, d, n);
        if (x == 1ULL || x == n - 1ULL) {
            return false;
        }
        for (int r = 1; r < s; ++r) {
            x = static_cast<u64>((static_cast<u128>(x) * x) % n);
            if (x == n - 1ULL) {
                return false;
            }
        }
        return true;
    };

    // Deterministic for 64-bit with this base set.
    for (const u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (witness(a)) {
            return false;
        }
    }
    return true;
}

u64 lagrange_from_0_to_m(const std::vector<u64>& y,
                         const int m,
                         const u64 x,
                         const std::vector<u64>& fac,
                         const std::vector<u64>& invfac,
                         const u64 p) {
    if (x <= static_cast<u64>(m)) {
        return y[static_cast<std::size_t>(x)];
    }

    std::vector<u64> pref(static_cast<std::size_t>(m + 2), 1ULL);
    std::vector<u64> suf(static_cast<std::size_t>(m + 2), 1ULL);

    for (int i = 0; i <= m; ++i) {
        pref[static_cast<std::size_t>(i + 1)] =
            static_cast<u64>((static_cast<u128>(pref[static_cast<std::size_t>(i)]) *
                              sub_mod(x, static_cast<u64>(i), p)) %
                             p);
    }
    for (int i = m; i >= 0; --i) {
        suf[static_cast<std::size_t>(i)] =
            static_cast<u64>((static_cast<u128>(suf[static_cast<std::size_t>(i + 1)]) *
                              sub_mod(x, static_cast<u64>(i), p)) %
                             p);
    }

    u64 out = 0ULL;
    for (int i = 0; i <= m; ++i) {
        u64 num =
            static_cast<u64>((static_cast<u128>(pref[static_cast<std::size_t>(i)]) *
                              suf[static_cast<std::size_t>(i + 1)]) %
                             p);
        u64 den =
            static_cast<u64>((static_cast<u128>(invfac[static_cast<std::size_t>(i)]) *
                              invfac[static_cast<std::size_t>(m - i)]) %
                             p);
        if (((m - i) & 1) != 0) {
            den = (den == 0ULL) ? 0ULL : (p - den);
        }
        u64 term =
            static_cast<u64>((static_cast<u128>(y[static_cast<std::size_t>(i)]) * num) % p);
        term = static_cast<u64>((static_cast<u128>(term) * den) % p);
        out = add_mod(out, term, p);
    }
    return out;
}

u64 power_sum_mod(const int exp, const u64 n, const u64 p, const std::vector<u64>& fac,
                  const std::vector<u64>& invfac) {
    const int m = exp + 1;
    std::vector<u64> y(static_cast<std::size_t>(m + 1), 0ULL);
    for (int i = 1; i <= m; ++i) {
        const u64 pw = mod_pow(static_cast<u64>(i), static_cast<u64>(exp), p);
        y[static_cast<std::size_t>(i)] = add_mod(y[static_cast<std::size_t>(i - 1)], pw, p);
    }
    return lagrange_from_0_to_m(y, m, n % p, fac, invfac, p);
}

u64 S_mod_prime(const int k, const u64 n, const u64 p) {
    const int max_m = k + 2;
    std::vector<u64> fac(static_cast<std::size_t>(max_m + 1), 1ULL);
    for (int i = 1; i <= max_m; ++i) {
        fac[static_cast<std::size_t>(i)] =
            static_cast<u64>((static_cast<u128>(fac[static_cast<std::size_t>(i - 1)]) * i) % p);
    }
    std::vector<u64> invfac(static_cast<std::size_t>(max_m + 1), 1ULL);
    invfac[static_cast<std::size_t>(max_m)] = mod_inv(fac[static_cast<std::size_t>(max_m)], p);
    for (int i = max_m; i >= 1; --i) {
        invfac[static_cast<std::size_t>(i - 1)] =
            static_cast<u64>((static_cast<u128>(invfac[static_cast<std::size_t>(i)] * i)) % p);
    }

    const u64 fk = power_sum_mod(k, n, p, fac, invfac);
    const u64 fk1 = power_sum_mod(k + 1, n, p, fac, invfac);
    const u64 n1 = (n + 1ULL) % p;
    const u64 part = static_cast<u64>((static_cast<u128>(n1) * fk) % p);
    return sub_mod(part, fk1, p);
}

u64 solve(const int k, const u64 n, const int p_lo, const int p_hi) {
    u64 total = 0ULL;
    for (int p = p_lo; p <= p_hi; ++p) {
        if (!is_prime_64(static_cast<u64>(p))) {
            continue;
        }
        total += S_mod_prime(k, n, static_cast<u64>(p));
    }
    return total;
}

u64 direct_S_small(const int k, const int n) {
    u128 out = 0;
    for (int i = 1; i <= n; ++i) {
        u128 fi = 0;
        for (int j = 1; j <= i; ++j) {
            u128 pw = 1;
            for (int t = 0; t < k; ++t) {
                pw *= static_cast<u128>(j);
            }
            fi += pw;
        }
        out += fi;
    }
    return static_cast<u64>(out);
}

bool run_checkpoints() {
    if (direct_S_small(4, 100) != 35'375'333'830ULL) {
        std::cerr << "Checkpoint failed: S_4(100)\n";
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

    std::cout << solve(options.k, options.n, options.p_lo, options.p_hi) << '\n';
    return 0;
}
