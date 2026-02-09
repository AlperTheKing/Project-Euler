#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i64 = std::int64_t;

struct Options {
    u64 n = 100000000000ULL;
    int k = 8;
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
    try {
        value = std::stoi(tail);
    } catch (...) {
        return false;
    }
    return true;
}

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }
    try {
        value = static_cast<u64>(std::stoull(tail));
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
        if (parse_u64_after_prefix(arg, "--n=", options.n) ||
            parse_int_after_prefix(arg, "--k=", options.k)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1ULL && options.k >= 1 && options.k <= 12;
}

u64 pow_u64(u64 base, int exp) {
    u64 r = 1;
    for (int i = 0; i < exp; ++i) {
        r *= base;
    }
    return r;
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0) {
        if (e & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        e >>= 1ULL;
    }
    return result;
}

i64 extended_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

u64 mod_inverse(u64 a, u64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    if (g != 1) {
        return 0ULL;
    }
    i64 r = x % static_cast<i64>(mod);
    if (r < 0) {
        r += static_cast<i64>(mod);
    }
    return static_cast<u64>(r);
}

std::vector<u64> roots_mod_power_10(int k) {
    const u64 m2 = 1ULL << k;
    const u64 m5 = pow_u64(5ULL, k);
    const u64 m = m2 * m5;

    std::vector<u64> roots2;
    for (u64 x = 0; x < m2; ++x) {
        if (mod_pow(x, 15ULL, m2) == (m2 - 1ULL) % m2) {
            roots2.push_back(x);
        }
    }

    std::vector<u64> roots5;
    for (u64 x = 0; x < m5; ++x) {
        if (mod_pow(x, 15ULL, m5) == (m5 - 1ULL) % m5) {
            roots5.push_back(x);
        }
    }

    const u64 inv_m2 = mod_inverse(m2 % m5, m5);

    std::vector<u64> residues;
    residues.reserve(roots2.size() * roots5.size());
    for (const u64 r2 : roots2) {
        for (const u64 r5 : roots5) {
            const u64 t = static_cast<u64>((static_cast<u128>((r5 + m5 - (r2 % m5)) % m5) * inv_m2) % m5);
            const u64 x = (r2 + m2 * t) % m;
            residues.push_back(x);
        }
    }

    std::sort(residues.begin(), residues.end());
    residues.erase(std::unique(residues.begin(), residues.end()), residues.end());
    return residues;
}

u128 sum_progression(u64 first, u64 step, u64 limit) {
    if (first > limit) {
        return 0;
    }
    const u64 cnt = (limit - first) / step + 1ULL;
    const u128 two_first = 2ULL * static_cast<u128>(first);
    const u128 last_term_part = (static_cast<u128>(cnt) - 1ULL) * static_cast<u128>(step);
    return static_cast<u128>(cnt) * (two_first + last_term_part) / 2ULL;
}

u128 f_value(u64 n, int k) {
    const u64 m = (1ULL << k) * pow_u64(5ULL, k);
    const std::vector<u64> residues = roots_mod_power_10(k);

    u128 total = 0;
    for (const u64 r : residues) {
        const u64 first = (r == 0ULL) ? m : r;
        total += sum_progression(first, m, n);
    }
    return total;
}

u128 brute_force(u64 n, u64 m) {
    u128 s = 0;
    for (u64 i = 1; i <= n; ++i) {
        if (mod_pow(i, 15ULL, m) == (m - 1ULL) % m) {
            s += i;
        }
    }
    return s;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }
    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool run_checkpoints() {
    if (f_value(100ULL, 2) != static_cast<u128>(295ULL)) {
        std::cerr << "Checkpoint failed: f(100,10^2)\n";
        return false;
    }
    if (f_value(1000ULL, 2) != static_cast<u128>(25450ULL)) {
        std::cerr << "Checkpoint failed: f(1000,10^2)\n";
        return false;
    }
    if (f_value(100000ULL, 4) != brute_force(100000ULL, 10000ULL)) {
        std::cerr << "Checkpoint failed: fast vs brute at m=10^4\n";
        return false;
    }
    const auto residues = roots_mod_power_10(8);
    if (residues.size() != 5U) {
        std::cerr << "Checkpoint failed: number of residues mod 10^8\n";
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

    std::cout << to_string_u128(f_value(options.n, options.k)) << '\n';
    return 0;
}
