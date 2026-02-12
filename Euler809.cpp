#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = static_cast<u64>((static_cast<u128>(r) * static_cast<u128>(a)) % static_cast<u128>(mod));
        }
        a = static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(a)) % static_cast<u128>(mod));
        e >>= 1ULL;
    }
    return r;
}

struct CycleData {
    u64 mu;
    u64 lam;
    std::vector<u64> values;
};

static CycleData build_s_cycle(u64 mod) {
    std::unordered_map<u64, u64> seen;
    seen.reserve(1 << 16);

    std::vector<u64> values;
    values.reserve(1 << 16);

    u64 x = 13 % mod;
    while (true) {
        auto it = seen.find(x);
        if (it != seen.end()) {
            const u64 mu = it->second;
            const u64 lam = static_cast<u64>(values.size()) - mu;
            return {mu, lam, values};
        }
        seen.emplace(x, static_cast<u64>(values.size()));
        values.push_back(x);
        const u64 y = pow_mod(2, x + 3, mod);
        x = (y + mod - 3) % mod;
    }
}

static u64 s_value_for_large_index(const CycleData& c, u64 index_mod_lam) {
    const u64 shift = c.mu % c.lam;
    const u64 k = (index_mod_lam + c.lam - shift) % c.lam;
    const u64 idx = c.mu + k;
    return c.values[idx];
}

static u64 inv_mod(u64 a, u64 mod) {
    std::int64_t t = 0;
    std::int64_t new_t = 1;
    std::int64_t r = static_cast<std::int64_t>(mod);
    std::int64_t new_r = static_cast<std::int64_t>(a % mod);

    while (new_r != 0) {
        const std::int64_t q = r / new_r;
        const std::int64_t tmp_t = t - q * new_t;
        t = new_t;
        new_t = tmp_t;
        const std::int64_t tmp_r = r - q * new_r;
        r = new_r;
        new_r = tmp_r;
    }

    assert(r == 1);
    if (t < 0) {
        t += static_cast<std::int64_t>(mod);
    }
    return static_cast<u64>(t);
}

int main() {
    const u64 F2 = 2;
    const u64 F3 = 3;
    const u64 F4 = 5;
    const u64 F5 = 13;
    const u64 F6 = (1ULL << (F5 + 3)) - 3;
    assert(F2 == 2);
    assert(F3 == 3);
    assert(F6 == 65533ULL);

    const u64 mod5 = 30'517'578'125ULL;
    const CycleData c5 = build_s_cycle(mod5);
    const CycleData c4 = build_s_cycle(c5.lam);
    const CycleData c3 = build_s_cycle(c4.lam);
    const CycleData c2 = build_s_cycle(c3.lam);

    assert(c5.mu == 25094 && c5.lam == 11824);
    assert(c4.mu == 0 && c4.lam == 17);
    assert(c3.mu == 0 && c3.lam == 3);
    assert(c2.mu == 0 && c2.lam == 1);

    u64 r5 = 65533 % mod5;
    u64 r4 = 65533 % c5.lam;
    u64 r3 = 65533 % c4.lam;
    u64 r2 = 65533 % c3.lam;

    auto advance_t_residues = [&]() {
        const u64 nr2 = s_value_for_large_index(c2, 0);
        const u64 nr3 = s_value_for_large_index(c3, r2);
        const u64 nr4 = s_value_for_large_index(c4, r3);
        const u64 nr5 = s_value_for_large_index(c5, r4);
        r2 = nr2;
        r3 = nr3;
        r4 = nr4;
        r5 = nr5;
    };

    advance_t_residues();
    advance_t_residues();
    const u64 stable5 = r5;
    advance_t_residues();
    assert(r5 == stable5);

    const u64 mod2 = 1ULL << 15;
    const u64 stable2 = mod2 - 3;

    const u64 inv = inv_mod(mod2 % mod5, mod5);
    const u64 delta = (stable5 + mod5 - stable2) % mod5;
    const u64 k = static_cast<u64>((static_cast<u128>(delta) * static_cast<u128>(inv)) % static_cast<u128>(mod5));
    const u64 answer = stable2 + mod2 * k;

    std::cout << answer << '\n';
    return 0;
}
