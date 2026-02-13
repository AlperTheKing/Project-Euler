#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u64 = std::uint64_t;

bool augment(int u, const std::vector<u64>& edges, std::vector<int>& match_right,
             std::vector<char>& seen) {
    u64 mask = edges[u];
    while (mask != 0) {
        const int v = __builtin_ctzll(mask);
        mask &= (mask - 1);
        if (seen[v]) {
            continue;
        }
        seen[v] = 1;
        if (match_right[v] == -1 || augment(match_right[v], edges, match_right, seen)) {
            match_right[v] = u;
            return true;
        }
    }
    return false;
}

bool has_perfect_matching(const std::vector<u64>& masks, int n) {
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return __builtin_popcountll(masks[a]) < __builtin_popcountll(masks[b]);
    });

    std::vector<u64> edges(n);
    for (int i = 0; i < n; ++i) {
        edges[i] = masks[order[i]];
        if (edges[i] == 0) {
            return false;
        }
    }

    std::vector<int> match_right(n, -1);
    for (int u = 0; u < n; ++u) {
        std::vector<char> seen(n, 0);
        if (!augment(u, edges, match_right, seen)) {
            return false;
        }
    }
    return true;
}

u64 lcm_upto(int n) {
    u64 l = 1;
    for (int i = 1; i <= n; ++i) {
        l = std::lcm(l, static_cast<u64>(i));
    }
    return l;
}

std::vector<int> prime_powers_upto(int n) {
    std::vector<int> factors;
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; p * p <= n; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        for (int q = p * p; q <= n; q += p) {
            is_prime[q] = false;
        }
    }
    for (int p = 2; p <= n; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        int pp = 1;
        while (pp <= n / p) {
            pp *= p;
        }
        factors.push_back(pp);
    }
    return factors;
}

std::vector<std::vector<std::vector<u64>>> build_mask_cache(int n) {
    std::vector<std::vector<std::vector<u64>>> cache(
        n + 1, std::vector<std::vector<u64>>(n + 1));
    for (int i = 1; i <= n; ++i) {
        for (int d = 1; d <= i; ++d) {
            if (i % d != 0) {
                continue;
            }
            cache[i][d].assign(d, 0);
            for (int r = 0; r < d; ++r) {
                u64 mask = 0;
                for (int o = 0; o < n; ++o) {
                    if ((r + o) % d == 0) {
                        mask |= (u64{1} << o);
                    }
                }
                cache[i][d][r] = mask;
            }
        }
    }
    return cache;
}

bool relaxed_feasible(int n, u64 modulus, u64 residue,
                      const std::vector<std::vector<std::vector<u64>>>& cache) {
    std::vector<u64> masks(n);
    for (int i = 1; i <= n; ++i) {
        const int d = static_cast<int>(std::gcd(modulus, static_cast<u64>(i)));
        const int r = static_cast<int>(residue % static_cast<u64>(d));
        masks[i - 1] = cache[i][d][r];
    }
    return has_perfect_matching(masks, n);
}

std::vector<u64> period_starts(int n) {
    const auto factors = prime_powers_upto(n);
    const auto cache = build_mask_cache(n);

    std::vector<u64> residues{0};
    u64 modulus = 1;

    for (int factor : factors) {
        const u64 next_modulus = modulus * static_cast<u64>(factor);
        std::vector<u64> next;
        next.reserve(residues.size() * static_cast<std::size_t>(factor));
        for (u64 c : residues) {
            for (int k = 0; k < factor; ++k) {
                const u64 candidate = c + static_cast<u64>(k) * modulus;
                if (relaxed_feasible(n, next_modulus, candidate, cache)) {
                    next.push_back(candidate);
                }
            }
        }
        residues.swap(next);
        modulus = next_modulus;
    }

    const u64 L = lcm_upto(n);
    assert(modulus == L);

    std::vector<u64> starts;
    starts.reserve(residues.size());
    for (u64 r : residues) {
        starts.push_back(r == 0 ? L : r);
    }
    std::sort(starts.begin(), starts.end());
    starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
    return starts;
}

bool is_divisible_range(u64 a, int n) {
    std::vector<u64> masks(n, 0);
    for (int i = 1; i <= n; ++i) {
        u64 mask = 0;
        for (int o = 0; o < n; ++o) {
            if ((a + static_cast<u64>(o)) % static_cast<u64>(i) == 0) {
                mask |= (u64{1} << o);
            }
        }
        masks[i - 1] = mask;
    }
    return has_perfect_matching(masks, n);
}

u64 nth_divisible_range_start(int n, u64 nth) {
    const u64 L = lcm_upto(n);
    const auto starts = period_starts(n);
    assert(!starts.empty());

    const u64 per = static_cast<u64>(starts.size());
    const u64 block = (nth - 1) / per;
    const u64 idx = (nth - 1) % per;
    return block * L + starts[static_cast<std::size_t>(idx)];
}

void validate() {
    std::vector<u64> first4;
    for (u64 a = 1; first4.size() < 4; ++a) {
        if (is_divisible_range(a, 4)) {
            first4.push_back(a);
        }
    }
    assert((first4 == std::vector<u64>{1, 2, 3, 6}));
    assert(nth_divisible_range_start(4, 4) == 6);

    assert(is_divisible_range(5, 36));
    assert(!is_divisible_range(6, 36));
}

}  // namespace

int main() {
    validate();
    std::cout << nth_divisible_range_start(36, 36) << '\n';
    return 0;
}
