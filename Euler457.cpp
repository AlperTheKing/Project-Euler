#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct Options {
    u32 l = 10'000'000U;
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
        out = static_cast<u32>(std::stoul(tail));
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
        if (parse_u32_after_prefix(arg, "--l=", options.l)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.l >= 2U;
}

u32 mod_pow(u64 base, u64 exp, const u32 mod) {
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

u32 tonelli_shanks(const u32 n, const u32 p) {
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

u64 lift_root(const u32 p, const u32 r) {
    const u32 deriv = static_cast<u32>((2ULL * r + p - 3ULL) % p);
    const u32 inv_deriv = mod_pow(deriv, p - 2U, p);

    const i128 fr = static_cast<i128>(r) * static_cast<i128>(r) - 3 * static_cast<i128>(r) - 1;
    const i64 q = static_cast<i64>(fr / static_cast<i128>(p));

    i64 neg_q = -(q % static_cast<i64>(p));
    neg_q %= static_cast<i64>(p);
    if (neg_q < 0) {
        neg_q += p;
    }

    const u32 t = static_cast<u32>((static_cast<u64>(neg_q) * inv_deriv) % p);
    return static_cast<u64>(r) + static_cast<u64>(t) * static_cast<u64>(p);
}

u128 solve(const u32 limit) {
    std::vector<bool> is_prime(static_cast<std::size_t>(limit) + 1U, true);
    if (limit >= 0U) {
        is_prime[0] = false;
    }
    if (limit >= 1U) {
        is_prime[1] = false;
    }

    for (u32 i = 2U; static_cast<u64>(i) * i <= limit; ++i) {
        if (!is_prime[i]) {
            continue;
        }
        for (u32 j = i * i; j <= limit; j += i) {
            is_prime[j] = false;
        }
    }

    std::array<bool, 13> residue{};
    residue.fill(false);
    residue[1] = true;
    residue[3] = true;
    residue[4] = true;
    residue[9] = true;
    residue[10] = true;
    residue[12] = true;

    u128 sum = 0;
    for (u32 p = 2U; p <= limit; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        if (p == 2U || p == 13U) {
            continue;
        }
        if (!residue[p % 13U]) {
            continue;
        }

        const u32 s = tonelli_shanks(13U % p, p);
        const u32 inv2 = (p + 1U) / 2U;

        const u32 r1 = static_cast<u32>((static_cast<u64>(3U + s) * inv2) % p);
        const u32 r2 = static_cast<u32>((static_cast<u64>(3U + p - s) * inv2) % p);

        const u64 n1 = lift_root(p, r1);
        const u64 n2 = lift_root(p, r2);
        sum += (n1 < n2 ? n1 : n2);
    }

    return sum;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const unsigned digit = static_cast<unsigned>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool run_checkpoints() {
    if (solve(10U) != 5U) {
        std::cerr << "Checkpoint failed: SR(10)\n";
        return false;
    }
    if (solve(100U) != 1'752U) {
        std::cerr << "Checkpoint failed: SR(100)\n";
        return false;
    }
    if (solve(1000U) != 6'728'355ULL) {
        std::cerr << "Checkpoint failed: SR(1000)\n";
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

    std::cout << to_string_u128(solve(options.l)) << '\n';
    return 0;
}
