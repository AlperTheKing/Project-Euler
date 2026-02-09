#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Options {
    u64 x = 1'000'000'000ULL;
    u64 y = 10'000'000ULL;
    u64 n = 1'000'000'000'000'000ULL;
    bool run_checkpoints = true;
};

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
        if (parse_u64_after_prefix(arg, "--x=", options.x)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--y=", options.y)) {
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }
    u64 out = 1ULL % mod;
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

u64 mod_inv_prime(const u64 a, const u64 p) { return mod_pow(a, p - 2ULL, p); }

u64 add_mod(const u64 a, const u64 b, const u64 p) {
    const u64 c = a + b;
    return (c >= p || c < a) ? (c % p) : c;
}

u64 sub_mod(const u64 a, const u64 b, const u64 p) {
    return (a >= b) ? (a - b) : (a + p - b);
}

u64 mul_mod(const u64 a, const u64 b, const u64 p) {
    return static_cast<u64>((static_cast<u128>(a) * b) % p);
}

std::vector<int> simple_primes_up_to(const int limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; static_cast<int64_t>(p) * p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int x = p * p; x <= limit; x += p) {
            is_prime[static_cast<std::size_t>(x)] = false;
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

std::vector<u64> primes_in_segment(const u64 lo, const u64 hi) {
    const int root = static_cast<int>(std::sqrt(static_cast<long double>(hi))) + 1;
    const std::vector<int> base_primes = simple_primes_up_to(root);

    std::vector<bool> is_prime(static_cast<std::size_t>(hi - lo + 1ULL), true);
    for (const int p : base_primes) {
        u64 start = (lo + static_cast<u64>(p) - 1ULL) / static_cast<u64>(p) * static_cast<u64>(p);
        if (start < static_cast<u64>(p) * static_cast<u64>(p)) {
            start = static_cast<u64>(p) * static_cast<u64>(p);
        }
        for (u64 x = start; x <= hi; x += static_cast<u64>(p)) {
            is_prime[static_cast<std::size_t>(x - lo)] = false;
        }
    }
    if (lo == 0ULL) {
        if (!is_prime.empty()) {
            is_prime[0] = false;
        }
        if (is_prime.size() > 1U) {
            is_prime[1] = false;
        }
    } else if (lo == 1ULL) {
        is_prime[0] = false;
    }

    std::vector<u64> out;
    out.reserve(static_cast<std::size_t>((hi - lo + 1ULL) / std::log(static_cast<long double>(hi))));
    for (u64 x = lo; x <= hi; ++x) {
        if (is_prime[static_cast<std::size_t>(x - lo)]) {
            out.push_back(x);
        }
    }
    return out;
}

u64 lucas_V_mod(const u64 k, const u64 p, const u64 P = 11ULL) {
    // Fast doubling on V_n(P,1):
    // V_{2m}   = V_m^2 - 2
    // V_{2m+1} = V_m V_{m+1} - P
    // V_{2m+2} = V_{m+1}^2 - 2
    if (k == 0ULL) {
        return 2ULL % p;
    }

    u64 vm = 2ULL % p;   // V_0
    u64 vmp1 = P % p;    // V_1

    int msb = 63 - __builtin_clzll(k);
    for (int bit = msb; bit >= 0; --bit) {
        const u64 v2m = sub_mod(mul_mod(vm, vm, p), 2ULL % p, p);
        const u64 v2m1 = sub_mod(mul_mod(vm, vmp1, p), P % p, p);
        const u64 v2m2 = sub_mod(mul_mod(vmp1, vmp1, p), 2ULL % p, p);

        if (((k >> bit) & 1ULL) == 0ULL) {
            vm = v2m;
            vmp1 = v2m1;
        } else {
            vm = v2m1;
            vmp1 = v2m2;
        }
    }
    return vm;
}

u64 a_n_mod_p(const u64 n, const u64 p) {
    if (n == 1ULL) {
        return 1ULL % p;
    }

    // b_n = 6 a_n + 5, b_{n+1}=b_n^2-2, b_1=11.
    // b_n = V_{2^{n-1}}(11,1) in F_p / F_{p^2}.
    // The index can be reduced modulo p-1 if 117 is a residue, otherwise p+1.
    const u64 leg = mod_pow(117ULL % p, (p - 1ULL) / 2ULL, p);
    const bool residue = (leg == 1ULL);
    const u64 order_mod = residue ? (p - 1ULL) : (p + 1ULL);

    const u64 idx = mod_pow(2ULL, n - 1ULL, order_mod);
    const u64 b = lucas_V_mod(idx, p, 11ULL);

    const u64 inv6 = mod_inv_prime(6ULL, p);
    const u64 a = mul_mod(sub_mod(b, 5ULL % p, p), inv6, p);
    return a;
}

u64 B(const u64 x, const u64 y, const u64 n) {
    const u64 lo = x;
    const u64 hi = x + y;
    const std::vector<u64> primes = primes_in_segment(lo, hi);

    u64 total = 0ULL;
    for (const u64 p : primes) {
        total += a_n_mod_p(n, p);
    }
    return total;
}

bool run_checkpoints() {
    if (a_n_mod_p(6ULL, 1'000'000'007ULL) != 203'064'689ULL) {
        std::cerr << "Checkpoint failed: a_6 mod 1e9+7\n";
        return false;
    }
    if (a_n_mod_p(100ULL, 1'000'000'007ULL) != 456'482'974ULL) {
        std::cerr << "Checkpoint failed: a_100 mod 1e9+7\n";
        return false;
    }
    if (B(1'000'000'000ULL, 1'000ULL, 1'000ULL) != 23'674'718'882ULL) {
        std::cerr << "Checkpoint failed: B(1e9,1e3,1e3)\n";
        return false;
    }
    if (B(1'000'000'000ULL, 1'000ULL, 1'000'000'000'000'000ULL) != 20'731'563'854ULL) {
        std::cerr << "Checkpoint failed: B(1e9,1e3,1e15)\n";
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

    std::cout << B(options.x, options.y, options.n) << '\n';
    return 0;
}
