#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

constexpr std::uint32_t kMod = 1'234'567'891U;

inline void add_mod(std::uint32_t& a, std::uint32_t b) {
    std::uint64_t s = static_cast<std::uint64_t>(a) + b;
    if (s >= kMod) {
        s -= kMod;
    }
    a = static_cast<std::uint32_t>(s);
}

inline void add_mul_mod(std::uint32_t& a, std::uint32_t x, std::uint32_t y) {
    const std::uint64_t prod = static_cast<std::uint64_t>(x) * static_cast<std::uint64_t>(y);
    std::uint64_t s = static_cast<std::uint64_t>(a) + (prod % kMod);
    if (s >= kMod) {
        s -= kMod;
    }
    a = static_cast<std::uint32_t>(s);
}

std::uint32_t solve(int N) {
    const int W = N + 1;
    auto idx = [W](int a, int b) -> std::size_t {
        return static_cast<std::size_t>(a) * static_cast<std::size_t>(W) + static_cast<std::size_t>(b);
    };

    std::vector<std::uint32_t> part(static_cast<std::size_t>(W) * static_cast<std::size_t>(W), 0);
    part[idx(0, 0)] = 1;
    for (int u = 1; u <= N; ++u) {
        part[idx(u, 0)] = 1;
        for (int v = 1; v <= N; ++v) {
            std::uint32_t val = part[idx(u - 1, v)];
            if (v >= u) {
                add_mod(val, part[idx(u, v - u)]);
            }
            part[idx(u, v)] = val;
        }
    }

    std::vector<std::uint32_t> g(static_cast<std::size_t>(W) * static_cast<std::size_t>(W), 0);
    for (int t = 0; t <= N; ++t) {
        for (int v = 0; v <= t; ++v) {
            g[idx(t, v)] = part[idx(t - v, v)];
        }
    }
    part.clear();
    part.shrink_to_fit();

    std::vector<std::uint32_t> agg(W, 0);
    std::vector<std::uint32_t> pref(W, 0);
    std::uint32_t contrib_ge2 = 0;

    for (int vA = 0; vA <= N; ++vA) {
        const int add_v = vA - 2;
        if (add_v >= 0) {
            for (int t = add_v; t <= N; ++t) {
                add_mod(agg[t], g[idx(t, add_v)]);
            }
        }

        std::uint32_t run = 0;
        for (int t = 0; t <= N; ++t) {
            add_mod(run, agg[t]);
            pref[t] = run;
        }

        for (int tA = vA; tA <= N; ++tA) {
            const std::uint32_t ga = g[idx(tA, vA)];
            if (ga == 0) {
                continue;
            }
            add_mul_mod(contrib_ge2, ga, pref[N - tA]);
        }
    }

    std::vector<std::uint32_t> pref_even(W, 0);
    std::vector<std::uint32_t> pref_odd(W, 0);
    std::uint32_t contrib_eq1 = 0;

    for (int vA = 1; vA <= N; ++vA) {
        const int vB = vA - 1;
        std::uint32_t even_sum = 0;
        std::uint32_t odd_sum = 0;

        for (int t = 0; t <= N; ++t) {
            const std::uint32_t gb = g[idx(t, vB)];
            if (t & 1) {
                add_mod(odd_sum, gb);
            } else {
                add_mod(even_sum, gb);
            }
            pref_even[t] = even_sum;
            pref_odd[t] = odd_sum;
        }

        for (int tA = vA; tA <= N; ++tA) {
            const std::uint32_t ga = g[idx(tA, vA)];
            if (ga == 0) {
                continue;
            }
            const int limit = N - tA;
            const std::uint32_t ways_b = (tA & 1) ? pref_odd[limit] : pref_even[limit];
            add_mul_mod(contrib_eq1, ga, ways_b);
        }
    }

    std::uint32_t ans = contrib_ge2;
    add_mod(ans, contrib_eq1);
    return ans;
}

void run_validations() {
    assert(solve(4) == 9U);
    assert(solve(10) == 486U);
    assert(solve(12) == 1309U);
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve(5000) << '\n';
    return 0;
}
