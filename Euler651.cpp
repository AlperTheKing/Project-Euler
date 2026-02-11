#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % kMod);
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) result = mod_mul(result, base);
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

std::vector<std::pair<u64, int>> factorize(u64 n) {
    std::vector<std::pair<u64, int>> factors;
    u64 x = n;
    for (u64 p = 2; p * p <= x; p += (p == 2 ? 1 : 2)) {
        if (x % p != 0) continue;
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        factors.emplace_back(p, e);
    }
    if (x > 1) factors.emplace_back(x, 1);
    return factors;
}

struct DivPhi {
    u64 d;
    u64 phi;
};

void gen_div_phi_rec(const std::vector<std::pair<u64, int>>& factors,
                     int idx,
                     u64 d,
                     u64 phi,
                     std::vector<DivPhi>& out) {
    if (idx == static_cast<int>(factors.size())) {
        out.push_back({d, phi});
        return;
    }

    const u64 p = factors[idx].first;
    const int emax = factors[idx].second;

    gen_div_phi_rec(factors, idx + 1, d, phi, out);

    u64 p_pow = 1;
    for (int e = 1; e <= emax; ++e) {
        p_pow *= p;
        u64 phi_factor = (p - 1);
        for (int i = 1; i < e; ++i) phi_factor *= p;
        gen_div_phi_rec(factors, idx + 1, d * p_pow, phi * phi_factor, out);
    }
}

std::vector<DivPhi> divisors_with_phi(u64 n) {
    std::vector<DivPhi> out;
    const auto factors = factorize(n);
    gen_div_phi_rec(factors, 0, 1, 1, out);
    return out;
}

struct CycleType {
    u64 length;
    u64 count;
};

struct ClassType {
    u64 multiplicity;
    std::vector<CycleType> cycles;
};

std::vector<ClassType> build_affine_classes(u64 n) {
    std::vector<ClassType> classes;

    if (n == 1) {
        classes.push_back({1, {{1, 1}}});
        return classes;
    }

    if (n == 2) {
        classes.push_back({1, {{1, 2}}});
        classes.push_back({1, {{2, 1}}});
        return classes;
    }

    const auto div_phi = divisors_with_phi(n);
    classes.reserve(div_phi.size() + 2);

    for (const auto& row : div_phi) {
        const u64 u = row.d;
        const u64 phi_u = row.phi;
        classes.push_back({phi_u, {{u, n / u}}});
    }

    if (n & 1ULL) {
        classes.push_back({n, {{1, 1}, {2, (n - 1) / 2}}});
    } else {
        classes.push_back({n / 2, {{1, 2}, {2, (n - 2) / 2}}});
        classes.push_back({n / 2, {{2, n / 2}}});
    }

    return classes;
}

u64 surjections_mod(u64 cycle_count, int colors, const std::array<std::array<u64, 41>, 41>& comb) {
    u64 result = 0;
    for (int j = 0; j <= colors; ++j) {
        const u64 ways_pick = comb[colors][j];
        const u64 ways_map = mod_pow(static_cast<u64>(colors - j), cycle_count);
        const u64 term = mod_mul(ways_pick, ways_map);
        if (j & 1) {
            result = (result + kMod - term) % kMod;
        } else {
            result += term;
            if (result >= kMod) result -= kMod;
        }
    }
    return result;
}

u64 cycle_count_product(const std::vector<CycleType>& x, const std::vector<CycleType>& y) {
    u128 total = 0;
    for (const auto& cx : x) {
        for (const auto& cy : y) {
            const u64 g = std::gcd(cx.length, cy.length);
            total += static_cast<u128>(cx.count) * static_cast<u128>(cy.count) * static_cast<u128>(g);
        }
    }
    return static_cast<u64>(total);
}

u64 solve_case(int colors, u64 a, u64 b, const std::array<std::array<u64, 41>, 41>& comb) {
    const auto classes_a = build_affine_classes(a);
    const auto classes_b = build_affine_classes(b);

    u64 group_a = 0;
    for (const auto& c : classes_a) group_a += c.multiplicity;
    u64 group_b = 0;
    for (const auto& c : classes_b) group_b += c.multiplicity;

    const u64 group_mod = mod_mul(group_a % kMod, group_b % kMod);
    const u64 inv_group = mod_pow(group_mod, kMod - 2);

    u64 burnside_sum = 0;
    std::unordered_map<u64, u64> cache;
    cache.reserve(classes_a.size() * classes_b.size() * 2);

    for (const auto& ca : classes_a) {
        for (const auto& cb : classes_b) {
            const u64 c = cycle_count_product(ca.cycles, cb.cycles);
            auto it = cache.find(c);
            u64 onto = 0;
            if (it == cache.end()) {
                onto = surjections_mod(c, colors, comb);
                cache.emplace(c, onto);
            } else {
                onto = it->second;
            }

            const u64 mult = mod_mul(ca.multiplicity % kMod, cb.multiplicity % kMod);
            burnside_sum += mod_mul(mult, onto);
            if (burnside_sum >= kMod) burnside_sum -= kMod;
        }
    }

    return mod_mul(burnside_sum, inv_group);
}

std::array<std::array<u64, 41>, 41> build_comb() {
    std::array<std::array<u64, 41>, 41> c{};
    c[0][0] = 1;
    for (int n = 1; n <= 40; ++n) {
        c[n][0] = 1;
        c[n][n] = 1;
        for (int k = 1; k < n; ++k) {
            c[n][k] = c[n - 1][k - 1] + c[n - 1][k];
        }
    }
    return c;
}

u64 solve() {
    const auto comb = build_comb();

    assert(solve_case(2, 2, 3, comb) == 11);
    assert(solve_case(3, 2, 3, comb) == 56);
    assert(solve_case(2, 3, 4, comb) == 156);
    assert(solve_case(8, 13, 21, comb) == 49718354);
    assert(solve_case(13, 144, 233, comb) == 907081451);

    std::vector<u64> fib(41, 0);
    fib[1] = 1;
    for (int i = 2; i <= 40; ++i) fib[i] = fib[i - 1] + fib[i - 2];

    u64 ans = 0;
    for (int i = 4; i <= 40; ++i) {
        ans += solve_case(i, fib[i - 1], fib[i], comb);
        ans %= kMod;
    }
    return ans;
}

}  // namespace

int main() {
    std::cout << solve() << "\n";
    return 0;
}
