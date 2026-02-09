#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using BigInt = boost::multiprecision::cpp_int;

constexpr i64 kMod = 1000000000LL;

struct Options {
    bool run_checkpoints = true;
};

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return true;
}

i64 mod_pow(i64 base, i64 exp, i64 mod) {
    i64 result = 1 % mod;
    i64 cur = base % mod;
    i64 e = exp;
    while (e > 0) {
        if (e & 1LL) {
            result = static_cast<i64>((__int128)result * cur % mod);
        }
        cur = static_cast<i64>((__int128)cur * cur % mod);
        e >>= 1LL;
    }
    return result;
}

i64 egcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

i64 mod_inv(i64 a, i64 mod) {
    i64 x = 0;
    i64 y = 0;
    const i64 g = egcd(a, mod, x, y);
    if (g != 1) {
        throw std::runtime_error("Inverse does not exist");
    }
    i64 res = x % mod;
    if (res < 0) {
        res += mod;
    }
    return res;
}

struct ModEntry {
    i64 mod;
    int pow2;
    int pow5;
    i64 period_mod;   // modulus needed for exponent (0 => not needed)
    i64 mod2_part;    // 2^pow2
    i64 mod5_part;    // 5^pow5
    i64 inv_mod2_in5; // inverse of mod2_part modulo mod5_part (if both parts exist)
};

std::vector<ModEntry> build_mod_chain() {
    std::vector<ModEntry> entries;

    auto pow5 = [](int k) -> i64 {
        i64 r = 1;
        for (int i = 0; i < k; ++i) {
            r *= 5;
        }
        return r;
    };

    // Main modulus 10^9 = 2^9 * 5^9, exponent period for odd part: 4*5^8.
    entries.push_back(ModEntry{kMod, 9, 9, 4 * pow5(8), 1LL << 9, pow5(9), 0});

    // Chain: 4*5^8, 4*5^7, ..., 4*5^1, 4.
    for (int e5 = 8; e5 >= 0; --e5) {
        const i64 mod = 4 * pow5(e5);
        const i64 period = (e5 == 0 ? 0 : 4 * pow5(e5 - 1));
        entries.push_back(ModEntry{mod, 2, e5, period, 4, pow5(e5), 0});
    }

    for (ModEntry& entry : entries) {
        if (entry.pow2 > 0 && entry.pow5 > 0) {
            entry.inv_mod2_in5 = mod_inv(entry.mod2_part % entry.mod5_part, entry.mod5_part);
        } else {
            entry.inv_mod2_in5 = 0;
        }
    }

    return entries;
}

i64 pow2_mod_composite(
    const ModEntry& entry,
    const std::vector<i64>& residues,
    const std::vector<ModEntry>& chain,
    const i64 exact_s,
    const bool exact_known) {
    if (entry.mod == 1) {
        return 0;
    }
    if (entry.pow5 == 0) {
        if (entry.pow2 == 0) {
            return 0;
        }
        if (exact_known && exact_s < entry.pow2) {
            return (1LL << exact_s) % entry.mod;
        }
        return 0;
    }

    i64 exp_mod = 0;
    if (entry.period_mod > 0) {
        int idx = -1;
        for (int i = 0; i < static_cast<int>(chain.size()); ++i) {
            if (chain[static_cast<std::size_t>(i)].mod == entry.period_mod) {
                idx = i;
                break;
            }
        }
        if (idx < 0) {
            throw std::runtime_error("Period modulus not in chain");
        }
        exp_mod = residues[static_cast<std::size_t>(idx)];
    }

    const i64 p5 = mod_pow(2, exp_mod, entry.mod5_part);

    if (entry.pow2 == 0) {
        return p5;
    }

    i64 p2 = 0;
    if (exact_known && exact_s < entry.pow2) {
        p2 = (1LL << exact_s) % entry.mod2_part;
    } else {
        p2 = 0;
    }

    // CRT for x ≡ p2 (mod 2^a), x ≡ p5 (mod 5^b).
    const i64 diff = (p5 - p2) % entry.mod5_part;
    const i64 adjusted = (diff + entry.mod5_part) % entry.mod5_part;
    const i64 t = static_cast<i64>((__int128)adjusted * entry.inv_mod2_in5 % entry.mod5_part);
    const i64 x = p2 + static_cast<i64>((__int128)entry.mod2_part * t % entry.mod);
    return x % entry.mod;
}

i64 compute_h_mod(i64 s0, int m) {
    const std::vector<ModEntry> chain = build_mod_chain();
    std::vector<i64> residues(chain.size(), 0);
    for (std::size_t i = 0; i < chain.size(); ++i) {
        residues[i] = s0 % chain[i].mod;
    }

    i64 exact_s = s0;
    bool exact_known = true;
    constexpr i64 kExactLimit = 1000;

    i64 sum_mod = 0;
    for (int i = 0; i <= m; ++i) {
        const i64 p_main = pow2_mod_composite(chain[0], residues, chain, exact_s, exact_known);
        const i64 term = static_cast<i64>((__int128)residues[0] * ((p_main + kMod - 1) % kMod) % kMod);
        sum_mod += term;
        if (sum_mod >= kMod) {
            sum_mod -= kMod;
        }

        if (i == m) {
            break;
        }

        std::vector<i64> next_residues(chain.size(), 0);
        for (std::size_t j = 0; j < chain.size(); ++j) {
            const i64 p = pow2_mod_composite(chain[j], residues, chain, exact_s, exact_known);
            next_residues[j] =
                static_cast<i64>((__int128)residues[j] * p % chain[j].mod);
        }
        residues.swap(next_residues);

        if (exact_known) {
            if (exact_s > 20) {
                exact_known = false;
            } else {
                i64 p2 = 1;
                for (int k = 0; k < exact_s; ++k) {
                    if (p2 > kExactLimit) {
                        break;
                    }
                    p2 <<= 1;
                }
                if (p2 > kExactLimit || exact_s > kExactLimit / std::max<i64>(1, p2)) {
                    exact_known = false;
                } else {
                    exact_s *= p2;
                    if (exact_s > kExactLimit) {
                        exact_known = false;
                    }
                }
            }
        }
    }

    return sum_mod;
}

i64 goodstein_small(const int n) {
    i64 g = n;
    int base = 2;
    i64 count = 0;
    while (g > 0) {
        ++count;

        std::vector<int> digits;
        i64 x = g;
        while (x > 0) {
            digits.push_back(static_cast<int>(x % base));
            x /= base;
        }
        if (digits.empty()) {
            digits.push_back(0);
        }

        i64 next = 0;
        for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
            next = next * static_cast<i64>(base + 1) + digits[static_cast<std::size_t>(i)];
        }
        g = next - 1;
        ++base;
    }
    return count;
}

i64 compute_g_mod(const int n, const std::array<i64, 8>& g_small) {
    if (n < 8) {
        return g_small[static_cast<std::size_t>(n)] % kMod;
    }
    const int low = n - 8;
    const i64 t = g_small[static_cast<std::size_t>(low)];
    const i64 s = t + 3;
    const int m = static_cast<int>(t + 2);
    const i64 h = compute_h_mod(s, m);
    return (t + h) % kMod;
}

i64 brute_h_mod_small(i64 s0, const int m) {
    BigInt s = s0;
    BigInt total = 0;
    const BigInt MOD = kMod;
    for (int i = 0; i <= m; ++i) {
        BigInt p = BigInt(1) << static_cast<unsigned>(s.convert_to<unsigned>());
        total += s * (p - 1);
        if (i < m) {
            s *= p;
        }
    }
    return static_cast<i64>((total % MOD).convert_to<i64>());
}

bool run_checkpoints() {
    if (goodstein_small(2) != 3) {
        std::cerr << "Checkpoint failed: G(2)\n";
        return false;
    }
    if (goodstein_small(4) != 21) {
        std::cerr << "Checkpoint failed: G(4)\n";
        return false;
    }
    if (goodstein_small(6) != 381) {
        std::cerr << "Checkpoint failed: G(6)\n";
        return false;
    }

    i64 sum_upto_7 = 0;
    for (int n = 1; n < 8; ++n) {
        sum_upto_7 += goodstein_small(n);
    }
    if (sum_upto_7 != 2517) {
        std::cerr << "Checkpoint failed: sum_{1<=n<8} G(n)\n";
        return false;
    }

    // Modular chain check on a still-manageable case with exact brute force.
    const i64 fast_h = compute_h_mod(2, 2);
    const i64 brute_h = brute_h_mod_small(2, 2);
    if (fast_h != brute_h) {
        std::cerr << "Checkpoint failed: H(s,m) modular engine\n";
        return false;
    }

    // Extra nontrivial modular checkpoint: H(3,2) can still be verified exactly.
    if (compute_h_mod(3, 2) != 722374141LL) {
        std::cerr << "Checkpoint failed: H(3,2)\n";
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

    std::array<i64, 8> g_small{};
    g_small[0] = 0;
    for (int n = 1; n < 8; ++n) {
        g_small[static_cast<std::size_t>(n)] = goodstein_small(n);
    }

    i64 answer = 0;
    for (int n = 1; n < 16; ++n) {
        answer += compute_g_mod(n, g_small);
        answer %= kMod;
    }

    std::cout << answer << '\n';
    return 0;
}
