#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    int limit = 20'000'000;
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
        std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_int_after_prefix(arg, "--limit=", options.limit)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.limit >= 3;
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

u32 mod_inverse(const u32 a, const u32 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    if (g != 1) {
        return 0U;
    }
    i64 v = x % static_cast<i64>(mod);
    if (v < 0) {
        v += static_cast<i64>(mod);
    }
    return static_cast<u32>(v);
}

std::vector<int> build_spf(const int limit) {
    std::vector<int> spf(static_cast<std::size_t>(limit + 1), 0);
    for (int i = 0; i <= limit; ++i) {
        spf[static_cast<std::size_t>(i)] = i;
    }
    for (int i = 2; i <= limit / i; ++i) {
        if (spf[static_cast<std::size_t>(i)] != i) {
            continue;
        }
        for (int j = i * i; j <= limit; j += i) {
            if (spf[static_cast<std::size_t>(j)] == j) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }
    return spf;
}

void roots_for_prime_power(const int prime, const int exponent, std::vector<u32>& roots, u32& modulus) {
    modulus = 1U;
    for (int i = 0; i < exponent; ++i) {
        modulus *= static_cast<u32>(prime);
    }

    roots.clear();
    if (prime == 2) {
        if (exponent == 1) {
            roots.push_back(1U);
        } else if (exponent == 2) {
            roots.push_back(1U);
            roots.push_back(3U);
        } else {
            const u32 half = 1U << static_cast<u32>(exponent - 1);
            roots.push_back(1U);
            roots.push_back(modulus - 1U);
            roots.push_back(1U + half);
            roots.push_back(modulus - 1U - half);
        }
    } else {
        roots.push_back(1U);
        roots.push_back(modulus - 1U);
    }
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
}

u32 compute_I(const int n, const std::vector<int>& spf) {
    int value = n;
    std::vector<std::pair<u32, std::vector<u32>>> blocks;
    blocks.reserve(8);

    while (value > 1) {
        const int p = spf[static_cast<std::size_t>(value)];
        int exp = 0;
        while (value % p == 0) {
            value /= p;
            ++exp;
        }
        u32 modulus = 0U;
        std::vector<u32> roots;
        roots_for_prime_power(p, exp, roots, modulus);
        blocks.push_back({modulus, std::move(roots)});
    }

    std::vector<u32> solutions;
    solutions.reserve(32);
    solutions.push_back(0U);
    u32 current_modulus = 1U;

    for (const auto& block : blocks) {
        const u32 component_mod = block.first;
        const std::vector<u32>& roots = block.second;
        const u32 inv = mod_inverse(current_modulus % component_mod, component_mod);

        std::vector<u32> next;
        next.reserve(solutions.size() * roots.size());
        for (const u32 x : solutions) {
            const u32 x_mod = x % component_mod;
            for (const u32 r : roots) {
                const u32 delta = (r >= x_mod) ? (r - x_mod) : (r + component_mod - x_mod);
                const u32 t = static_cast<u32>((static_cast<u64>(delta) * inv) % component_mod);
                const u32 candidate = x + current_modulus * t;
                next.push_back(candidate);
            }
        }
        solutions.swap(next);
        current_modulus *= component_mod;
    }

    u32 best = 1U;
    const u32 target = static_cast<u32>(n - 1);
    for (const u32 x : solutions) {
        if (x < target && x > best) {
            best = x;
        }
    }
    return best;
}

u64 solve(const int limit) {
    const std::vector<int> spf = build_spf(limit);
    u64 sum = 0ULL;
    for (int n = 3; n <= limit; ++n) {
        sum += compute_I(n, spf);
    }
    return sum;
}

u32 brute_I(const int n) {
    u32 best = 1U;
    for (int m = 1; m < n - 1; ++m) {
        if (std::gcd(m, n) != 1) {
            continue;
        }
        if ((static_cast<u64>(m) * static_cast<u64>(m)) % static_cast<u64>(n) == 1ULL) {
            best = static_cast<u32>(m);
        }
    }
    return best;
}

bool run_checkpoints() {
    const std::vector<int> spf = build_spf(2'000);
    if (compute_I(7, spf) != 1U) {
        std::cerr << "Checkpoint failed: I(7)=1" << '\n';
        return false;
    }
    if (compute_I(100, spf) != 51U) {
        std::cerr << "Checkpoint failed: I(100)=51" << '\n';
        return false;
    }
    for (int n = 3; n <= 600; ++n) {
        if (compute_I(n, spf) != brute_I(n)) {
            std::cerr << "Checkpoint failed: brute-force cross-check at n=" << n << '\n';
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
    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }
    std::cout << solve(options.limit) << '\n';
    return 0;
}
