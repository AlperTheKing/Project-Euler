#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 123'456'789ULL;

u64 next_value(u64 x, u64 mod) {
    const u64 y = (x + mod - 1) % mod;
    const u64 y2 = static_cast<u64>((u128)y * y % mod);
    const u64 y3 = static_cast<u64>((u128)y2 * y % mod);
    return (y3 + 2) % mod;
}

std::pair<u64, u64> cycle_mu_lambda(u64 mod) {
    u64 tortoise = next_value(1, mod);
    u64 hare = next_value(next_value(1, mod), mod);

    while (tortoise != hare) {
        tortoise = next_value(tortoise, mod);
        hare = next_value(next_value(hare, mod), mod);
    }

    u64 mu = 0;
    tortoise = 1;
    while (tortoise != hare) {
        tortoise = next_value(tortoise, mod);
        hare = next_value(hare, mod);
        ++mu;
    }

    u64 lambda = 1;
    hare = next_value(tortoise, mod);
    while (tortoise != hare) {
        hare = next_value(hare, mod);
        ++lambda;
    }

    return {mu, lambda};
}

struct Group {
    u32 l;
    u32 r;
    u32 q;
};

u64 solve(u32 N) {
    const auto [mu0, lambda] = cycle_mu_lambda(kMod);
    const u64 cycle_start = mu0 + 1;

    const u32 need = static_cast<u32>(cycle_start + lambda);
    std::vector<u32> s_mod(need + 1, 0);
    s_mod[1] = 1;
    for (u32 i = 2; i <= need; ++i) {
        s_mod[i] = static_cast<u32>(next_value(s_mod[i - 1], kMod));
    }

    auto value_by_index = [&](u64 k) -> u32 {
        if (k < cycle_start) {
            return s_mod[static_cast<u32>(k)];
        }
        const u64 idx = cycle_start + (k - cycle_start) % lambda;
        return s_mod[static_cast<u32>(idx)];
    };

    std::vector<Group> groups;
    groups.reserve(2 * static_cast<std::size_t>(std::sqrt(static_cast<long double>(N))) + 8);

    for (u32 l = 1; l <= N;) {
        const u32 q = N / l;
        const u32 r = N / q;
        groups.push_back({l, r, q});
        l = r + 1;
    }

    std::vector<u32> q_values;
    q_values.reserve(groups.size());
    for (const auto& g : groups) {
        q_values.push_back(g.q);
    }
    std::sort(q_values.begin(), q_values.end());
    q_values.erase(std::unique(q_values.begin(), q_values.end()), q_values.end());

    std::unordered_map<u32, u32> coprime_pair_count_mod;
    coprime_pair_count_mod.reserve(q_values.size() * 2);

    std::vector<u32> phi(static_cast<std::size_t>(N) + 1);
    for (u32 i = 0; i <= N; ++i) {
        phi[i] = i;
    }

    for (u32 p = 2; p <= N; ++p) {
        if (phi[p] != p) {
            continue;
        }
        for (u32 j = p; j <= N; j += p) {
            phi[j] -= phi[j] / p;
        }
    }

    i64 prefix_phi_mod = 0;
    std::size_t q_ptr = 0;
    for (u32 i = 1; i <= N; ++i) {
        prefix_phi_mod += phi[i];
        prefix_phi_mod %= static_cast<i64>(kMod);

        while (q_ptr < q_values.size() && q_values[q_ptr] == i) {
            i64 v = (2 * prefix_phi_mod - 1) % static_cast<i64>(kMod);
            if (v < 0) {
                v += static_cast<i64>(kMod);
            }
            coprime_pair_count_mod[q_values[q_ptr]] = static_cast<u32>(v);
            ++q_ptr;
        }
    }

    phi.clear();
    phi.shrink_to_fit();

    i64 answer = 0;
    i64 pref_u_mod = 0;
    i64 prev_pref_u_mod = 0;

    u32 d = 1;
    u64 s_d_mod_lambda = 1 % lambda;

    constexpr u64 small_index[5] = {0, 1, 2, 3, 10};

    for (const auto& g : groups) {
        while (d <= g.r) {
            u32 u_d;
            if (d <= 4) {
                u_d = value_by_index(small_index[d]);
            } else {
                const u64 start_mod = cycle_start % lambda;
                const u64 delta = (s_d_mod_lambda + lambda - start_mod) % lambda;
                const u64 idx = cycle_start + delta;
                u_d = s_mod[static_cast<u32>(idx)];
            }

            pref_u_mod += u_d;
            pref_u_mod %= static_cast<i64>(kMod);

            s_d_mod_lambda = next_value(s_d_mod_lambda, lambda);
            ++d;
        }

        i64 segment = pref_u_mod - prev_pref_u_mod;
        segment %= static_cast<i64>(kMod);
        if (segment < 0) {
            segment += static_cast<i64>(kMod);
        }

        const u32 cmod = coprime_pair_count_mod[g.q];
        answer += static_cast<i64>((u128)segment * cmod % kMod);
        answer %= static_cast<i64>(kMod);

        prev_pref_u_mod = pref_u_mod;
    }

    if (answer < 0) {
        answer += static_cast<i64>(kMod);
    }
    return static_cast<u64>(answer);
}

void validate() {
    assert(solve(3) == 12);
    assert(solve(4) == 24'881'925ULL);
    assert(solve(100) == 14'416'749ULL);
}

}  // namespace

int main() {
    validate();
    std::cout << solve(100'000'000) << '\n';
    return 0;
}
