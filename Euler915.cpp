#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <unordered_map>

namespace {

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 123'456'789ULL;

u64 next_value(u64 x, u64 mod) {
    const u64 y = (x == 0 ? mod - 1 : x - 1);
    const u64 y2 = static_cast<u64>((u128)y * y % mod);
    const u64 y3 = static_cast<u64>((u128)y2 * y % mod);
    u64 z = y3 + 2;
    if (z >= mod) {
        z -= mod;
    }
    return z;
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

class SumPhi {
public:
    explicit SumPhi(u32 n) {
        constexpr u32 kBase = 1'000'000;
        limit_ = (n < kBase ? n : kBase);

        std::vector<u32> lp(static_cast<std::size_t>(limit_) + 1, 0);
        std::vector<u32> phi(static_cast<std::size_t>(limit_) + 1, 0);
        phi[1] = 1;

        std::vector<u32> primes;
        primes.reserve(static_cast<std::size_t>(limit_ / 10));
        for (u32 i = 2; i <= limit_; ++i) {
            if (lp[i] == 0) {
                lp[i] = i;
                primes.push_back(i);
                phi[i] = i - 1;
            }
            for (u32 p : primes) {
                u64 v = static_cast<u64>(i) * p;
                if (v > limit_ || p > lp[i]) {
                    break;
                }
                lp[static_cast<u32>(v)] = p;
                if (i % p == 0) {
                    phi[static_cast<u32>(v)] = phi[i] * p;
                    break;
                }
                phi[static_cast<u32>(v)] = phi[i] * (p - 1);
            }
        }

        prefix_.assign(static_cast<std::size_t>(limit_) + 1, 0);
        for (u32 i = 1; i <= limit_; ++i) {
            prefix_[i] = prefix_[i - 1] + phi[i];
        }

        memo_.reserve(1 << 17);
    }

    u64 get(u64 n) {
        if (n <= limit_) {
            return prefix_[static_cast<std::size_t>(n)];
        }
        const auto it = memo_.find(n);
        if (it != memo_.end()) {
            return it->second;
        }

        u64 ans = static_cast<u64>((static_cast<u128>(n) * (n + 1)) / 2);
        for (u64 l = 2; l <= n;) {
            const u64 q = n / l;
            const u64 r = n / q;
            const u64 cnt = r - l + 1;
            ans -= cnt * get(q);
            l = r + 1;
        }
        memo_[n] = ans;
        return ans;
    }

private:
    u32 limit_ = 0;
    std::vector<u64> prefix_;
    std::unordered_map<u64, u64> memo_;
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

    SumPhi sum_phi(N);

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

        const u64 sp = sum_phi.get(g.q) % kMod;
        const u32 cmod = static_cast<u32>((2ULL * sp + kMod - 1ULL) % kMod);
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
