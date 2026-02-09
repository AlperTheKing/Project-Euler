#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u32 kMod = 1'000'000'007U;

struct Options {
    u32 n = 10'000'000U;
    bool run_checkpoints = true;
};

bool parse_u32_after_prefix(const std::string& arg, const std::string& prefix, u32& out) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        const unsigned long long parsed = std::stoull(tail);
        if (parsed > static_cast<unsigned long long>(std::numeric_limits<u32>::max())) {
            return false;
        }
        out = static_cast<u32>(parsed);
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
        if (parse_u32_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1U;
}

u32 mod_pow(u64 base, u64 exp, const u32 mod = kMod) {
    u64 result = 1ULL;
    base %= mod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1ULL;
    }
    return static_cast<u32>(result);
}

u32 tonelli_shanks(u32 n, u32 p) {
    if (p == 2U) {
        return n & 1U;
    }
    if (n == 0U) {
        return 0U;
    }
    if (mod_pow(n, (p - 1U) / 2U, p) != 1U) {
        return 0U;
    }
    if ((p & 3U) == 3U) {
        return mod_pow(n, (p + 1U) / 4U, p);
    }

    u32 q = p - 1U;
    u32 s = 0U;
    while ((q & 1U) == 0U) {
        q >>= 1U;
        ++s;
    }

    u32 z = 2U;
    while (mod_pow(z, (p - 1U) / 2U, p) != p - 1U) {
        ++z;
    }

    u64 m = s;
    u64 c = mod_pow(z, q, p);
    u64 t = mod_pow(n, q, p);
    u64 r = mod_pow(n, (q + 1U) / 2U, p);

    while (t != 1U) {
        u64 tt = t;
        u64 i = 0U;
        while (tt != 1U && i < m) {
            tt = (tt * tt) % p;
            ++i;
        }

        const u64 shift = m - i - 1U;
        const u64 b = mod_pow(static_cast<u32>(c), 1ULL << shift, p);
        r = (r * b) % p;
        const u64 b2 = (b * b) % p;
        t = (t * b2) % p;
        c = b2;
        m = i;
    }

    return static_cast<u32>(r);
}

u64 brute_f_1024_exact() {
    auto r_exact = [](u64 n) -> u64 {
        u64 m = n;
        u64 q = 1ULL;
        for (u64 p = 2ULL; p * p <= m; ++p) {
            if (m % p != 0ULL) {
                continue;
            }
            u64 pe = 1ULL;
            while (m % p == 0ULL) {
                m /= p;
                pe *= p;
            }
            q *= (1ULL + pe);
        }
        if (m > 1ULL) {
            q *= (1ULL + m);
        }
        return q - n;
    };

    u64 sum = 0ULL;
    for (u64 n = 1ULL; n <= 1024ULL; ++n) {
        const u64 v = n * n * n * n + 4ULL;
        sum += r_exact(v);
    }
    return sum;
}

u32 solve_mod(const u32 n_limit) {
    const u32 max_m = n_limit + 1U;
    std::vector<u64> remaining(static_cast<std::size_t>(max_m) + 1U, 0ULL);
    std::vector<u32> q_value(static_cast<std::size_t>(max_m) + 1U, 1U);

    for (u32 m = 0U; m <= max_m; ++m) {
        remaining[m] = static_cast<u64>(m) * static_cast<u64>(m) + 1ULL;
    }

    std::vector<u32> spf(static_cast<std::size_t>(max_m) + 1U, 0U);
    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(max_m / 10U));

    for (u32 i = 2U; i <= max_m; ++i) {
        if (spf[i] == 0U) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (const u32 p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > max_m || p > spf[i]) {
                break;
            }
            spf[static_cast<std::size_t>(v)] = p;
        }
    }

    const auto apply_prime_to_index = [&](const u32 p, const u32 idx) {
        u64 value = remaining[idx];
        if (value % p != 0ULL) {
            return;
        }
        u64 p_power_mod = 1ULL;
        do {
            value /= p;
            p_power_mod = (p_power_mod * p) % kMod;
        } while (value % p == 0ULL);
        remaining[idx] = value;
        const u32 term = static_cast<u32>((1ULL + p_power_mod) % kMod);
        q_value[idx] = static_cast<u32>((static_cast<u64>(q_value[idx]) * term) % kMod);
    };

    for (const u32 p : primes) {
        if (p == 2U) {
            for (u32 idx = 1U; idx <= max_m; idx += 2U) {
                apply_prime_to_index(2U, idx);
            }
            continue;
        }
        if ((p & 3U) != 1U) {
            continue;
        }

        const u32 root = tonelli_shanks(p - 1U, p);
        if (root == 0U) {
            continue;
        }
        const u32 root2 = (p - root) % p;

        for (u32 idx = root; idx <= max_m; idx += p) {
            apply_prime_to_index(p, idx);
        }
        if (root2 != root) {
            for (u32 idx = root2; idx <= max_m; idx += p) {
                apply_prime_to_index(p, idx);
            }
        }
    }

    for (u32 m = 0U; m <= max_m; ++m) {
        const u64 leftover = remaining[m];
        if (leftover > 1ULL) {
            const u32 term = static_cast<u32>((1ULL + (leftover % kMod)) % kMod);
            q_value[m] = static_cast<u32>((static_cast<u64>(q_value[m]) * term) % kMod);
        }
    }

    const u32 inv9 = mod_pow(9U, static_cast<u64>(kMod) - 2ULL);
    const u32 even_adjust = static_cast<u32>((5ULL * inv9) % kMod);

    u64 answer = 0ULL;
    for (u32 n = 1U; n <= n_limit; ++n) {
        const u32 q1 = q_value[n - 1U];
        const u32 q2 = q_value[n + 1U];
        u64 q_total = (static_cast<u64>(q1) * q2) % kMod;
        if ((n & 1U) == 0U) {
            q_total = (q_total * even_adjust) % kMod;
        }

        const u64 n_mod = n % kMod;
        const u64 n2 = (n_mod * n_mod) % kMod;
        const u64 n4_plus_4 = (n2 * n2 + 4ULL) % kMod;

        const u64 r_mod = (q_total + kMod - n4_plus_4) % kMod;
        answer += r_mod;
        answer %= kMod;
    }

    return static_cast<u32>(answer);
}

bool run_checkpoints() {
    if (brute_f_1024_exact() != 77'532'377'300'600ULL) {
        std::cerr << "Checkpoint failed: exact F(1024)\n";
        return false;
    }
    if (solve_mod(1024U) != 376'757'876U) {
        std::cerr << "Checkpoint failed: F(1024) mod 1e9+7\n";
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

    std::cout << solve_mod(options.n) << '\n';
    return 0;
}
