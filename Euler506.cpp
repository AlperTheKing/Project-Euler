#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 123'454'321ULL;
constexpr u64 kBlockBase = 1'000'000ULL;  // 10^6
constexpr std::array<int, 6> kDigits = {1, 2, 3, 4, 3, 2};

struct DigitStream {
    int pos = 0;

    int next() {
        const int d = kDigits[static_cast<std::size_t>(pos)];
        pos = (pos + 1) % static_cast<int>(kDigits.size());
        return d;
    }
};

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        e >>= 1ULL;
    }
    return result;
}

std::int64_t extended_gcd(const std::int64_t a, const std::int64_t b, std::int64_t& x,
                          std::int64_t& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    std::int64_t x1 = 0;
    std::int64_t y1 = 0;
    const std::int64_t g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

u64 mod_inverse(const u64 a, const u64 mod) {
    std::int64_t x = 0;
    std::int64_t y = 0;
    const std::int64_t g = extended_gcd(static_cast<std::int64_t>(a),
                                        static_cast<std::int64_t>(mod), x, y);
    if (g != 1) {
        return 0ULL;
    }
    const std::int64_t m = static_cast<std::int64_t>(mod);
    std::int64_t res = x % m;
    if (res < 0) {
        res += m;
    }
    return static_cast<u64>(res);
}

std::vector<u128> direct_terms_exact(const int count) {
    DigitStream stream;
    std::vector<u128> out(static_cast<std::size_t>(count + 1), 0);
    for (int n = 1; n <= count; ++n) {
        int sum = 0;
        u128 value = 0;
        while (sum < n) {
            const int d = stream.next();
            sum += d;
            value = value * 10 + static_cast<u128>(d);
        }
        if (sum != n) {
            return {};
        }
        out[static_cast<std::size_t>(n)] = value;
    }
    return out;
}

u64 direct_sum_mod(const int n_max) {
    DigitStream stream;
    u64 total = 0ULL;
    for (int n = 1; n <= n_max; ++n) {
        int sum = 0;
        u64 value_mod = 0ULL;
        while (sum < n) {
            const int d = stream.next();
            sum += d;
            value_mod = (value_mod * 10ULL + static_cast<u64>(d)) % kMod;
        }
        if (sum != n) {
            return 0ULL;
        }
        total += value_mod;
        total %= kMod;
    }
    return total;
}

u64 direct_sum_exact_small(const int n_max) {
    DigitStream stream;
    u64 total = 0ULL;
    for (int n = 1; n <= n_max; ++n) {
        int sum = 0;
        u64 value = 0ULL;
        while (sum < n) {
            const int d = stream.next();
            sum += d;
            value = value * 10ULL + static_cast<u64>(d);
        }
        if (sum != n) {
            return 0ULL;
        }
        total += value;
    }
    return total;
}

u64 geometric_sum(const u64 ratio, const u64 terms, const u64 mod, const u64 inv_ratio_minus_one) {
    if (terms == 0ULL) {
        return 0ULL;
    }
    const u64 p = mod_pow(ratio, terms, mod);
    const u64 numerator = (p + mod - 1ULL) % mod;
    return static_cast<u64>((static_cast<u128>(numerator) * inv_ratio_minus_one) % mod);
}

u64 fast_sum_mod(const u64 n_max) {
    const std::vector<u128> terms = direct_terms_exact(45);
    if (terms.empty()) {
        return 0ULL;
    }

    std::array<u64, 15> base{};
    std::array<u64, 15> add{};
    for (int r = 1; r <= 15; ++r) {
        const u128 v0 = terms[static_cast<std::size_t>(r)];
        const u128 v1 = terms[static_cast<std::size_t>(r + 15)];
        const u128 v2 = terms[static_cast<std::size_t>(r + 30)];
        const u128 c = v1 - v0 * static_cast<u128>(kBlockBase);
        if (v2 != v1 * static_cast<u128>(kBlockBase) + c) {
            return 0ULL;
        }
        base[static_cast<std::size_t>(r - 1)] = static_cast<u64>(v0);
        add[static_cast<std::size_t>(r - 1)] = static_cast<u64>(c);
    }

    const u64 inv_b_minus_1 = mod_inverse(kBlockBase - 1ULL, kMod);
    if (inv_b_minus_1 == 0ULL) {
        return 0ULL;
    }

    u64 total = 0ULL;
    for (int r = 1; r <= 15; ++r) {
        if (n_max < static_cast<u64>(r)) {
            break;
        }
        const u64 count = (n_max - static_cast<u64>(r)) / 15ULL + 1ULL;
        const u64 g = geometric_sum(kBlockBase, count, kMod, inv_b_minus_1);
        const u64 count_mod = count % kMod;
        const u64 g_minus_count = (g + kMod - count_mod) % kMod;
        const u64 t = static_cast<u64>((static_cast<u128>(g_minus_count) * inv_b_minus_1) % kMod);

        const u64 part_base =
            static_cast<u64>((static_cast<u128>(base[static_cast<std::size_t>(r - 1)] % kMod) * g) %
                             kMod);
        const u64 part_add =
            static_cast<u64>((static_cast<u128>(add[static_cast<std::size_t>(r - 1)] % kMod) * t) %
                             kMod);
        total += part_base;
        total %= kMod;
        total += part_add;
        total %= kMod;
    }
    return total;
}

bool run_checkpoints() {
    if (direct_sum_exact_small(11) != 36'120ULL) {
        std::cerr << "Checkpoint failed: S(11)\n";
        return false;
    }
    if (direct_sum_mod(1'000) != 18'232'686ULL) {
        std::cerr << "Checkpoint failed: S(1000) mod M\n";
        return false;
    }
    if (fast_sum_mod(5'000ULL) != direct_sum_mod(5'000)) {
        std::cerr << "Checkpoint failed: fast/direct consistency at 5000\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 100'000'000'000'000ULL;
    std::cout << fast_sum_mod(n) << '\n';
    return 0;
}
