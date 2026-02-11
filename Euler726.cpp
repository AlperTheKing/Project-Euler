#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using u64 = std::uint64_t;
using boost::multiprecision::cpp_int;

constexpr u64 kMod = 1'000'000'033ULL;

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = (result * base) % kMod;
        }
        base = (base * base) % kMod;
        exp >>= 1ULL;
    }
    return result;
}

u64 mod_inv(const u64 x) {
    return mod_pow(x, kMod - 2ULL);
}

std::vector<u64> compute_f_mod(const int limit) {
    std::vector<u64> f(static_cast<std::size_t>(limit + 1), 0ULL);
    f[0] = 1ULL;

    u64 prefix_a = 1ULL;
    u64 power_part = 1ULL;
    u64 factorial_tri = 1ULL;
    u64 pow2 = 1ULL;
    u64 tri = 0ULL;

    for (int n = 1; n <= limit; ++n) {
        pow2 = (pow2 * 2ULL) % kMod;
        const u64 numer = (pow2 + kMod - 1ULL) % kMod;
        const u64 denom_inv = mod_inv(static_cast<u64>(2 * n - 1));
        const u64 a_n = (numer * denom_inv) % kMod;

        prefix_a = (prefix_a * a_n) % kMod;
        power_part = (power_part * prefix_a) % kMod;

        for (int i = 0; i < n; ++i) {
            ++tri;
            factorial_tri = (factorial_tri * tri) % kMod;
        }

        f[static_cast<std::size_t>(n)] = (factorial_tri * power_part) % kMod;
    }

    return f;
}

struct TinySolver {
    int n = 0;
    int node_count = 0;
    std::vector<int> row_of;
    std::vector<std::vector<int>> parents;
    std::vector<int> topo;
    std::unordered_map<u64, cpp_int> memo;

    explicit TinySolver(const int n_layers) : n(n_layers) {
        node_count = n * (n + 1) / 2;
        row_of.assign(static_cast<std::size_t>(node_count), 0);
        parents.assign(static_cast<std::size_t>(node_count), {});
        topo.reserve(static_cast<std::size_t>(node_count));

        std::vector<std::vector<int>> id(static_cast<std::size_t>(n + 1));
        int idx = 0;
        for (int r = 1; r <= n; ++r) {
            id[static_cast<std::size_t>(r)].resize(static_cast<std::size_t>(r));
            for (int c = 1; c <= r; ++c) {
                id[static_cast<std::size_t>(r)][static_cast<std::size_t>(c - 1)] = idx;
                row_of[static_cast<std::size_t>(idx)] = r;
                topo.push_back(idx);
                ++idx;
            }
        }

        for (int r = 1; r <= n; ++r) {
            for (int c = 1; c <= r; ++c) {
                const int v = id[static_cast<std::size_t>(r)][static_cast<std::size_t>(c - 1)];
                auto& p = parents[static_cast<std::size_t>(v)];
                if (r > 1 && c > 1) {
                    p.push_back(id[static_cast<std::size_t>(r - 1)][static_cast<std::size_t>(c - 2)]);
                }
                if (r > 1 && c < r) {
                    p.push_back(id[static_cast<std::size_t>(r - 1)][static_cast<std::size_t>(c - 1)]);
                }
            }
        }
    }

    cpp_int solve_mask(const u64 mask) {
        if (mask == 0ULL) {
            return 1;
        }
        const auto it = memo.find(mask);
        if (it != memo.end()) {
            return it->second;
        }

        std::vector<int> maxima;
        maxima.reserve(static_cast<std::size_t>(n));
        for (int v = 0; v < node_count; ++v) {
            if (((mask >> v) & 1ULL) == 0ULL) {
                continue;
            }
            bool has_parent = false;
            for (const int p : parents[static_cast<std::size_t>(v)]) {
                if (((mask >> p) & 1ULL) != 0ULL) {
                    has_parent = true;
                    break;
                }
            }
            if (!has_parent) {
                maxima.push_back(v);
            }
        }

        cpp_int result = 0;
        for (const int m : maxima) {
            std::vector<cpp_int> ways(static_cast<std::size_t>(node_count), 0);
            ways[static_cast<std::size_t>(m)] = 1;

            for (const int v : topo) {
                if (row_of[static_cast<std::size_t>(v)] <= row_of[static_cast<std::size_t>(m)]) {
                    continue;
                }
                if (((mask >> v) & 1ULL) == 0ULL) {
                    continue;
                }
                cpp_int count = 0;
                for (const int p : parents[static_cast<std::size_t>(v)]) {
                    if (((mask >> p) & 1ULL) != 0ULL) {
                        count += ways[static_cast<std::size_t>(p)];
                    }
                }
                ways[static_cast<std::size_t>(v)] = count;
            }

            cpp_int weight = 0;
            for (int v = 0; v < node_count; ++v) {
                if (((mask >> v) & 1ULL) != 0ULL) {
                    weight += ways[static_cast<std::size_t>(v)];
                }
            }

            result += weight * solve_mask(mask & ~(1ULL << m));
        }

        memo.emplace(mask, result);
        return result;
    }

    cpp_int solve_full() {
        const u64 full_mask = (node_count == 64) ? ~0ULL : ((1ULL << node_count) - 1ULL);
        return solve_mask(full_mask);
    }
};

u64 mod_from_cpp_int(const cpp_int& value) {
    cpp_int x = value % kMod;
    if (x < 0) {
        x += kMod;
    }
    return x.convert_to<u64>();
}

}  // namespace

int main() {
    {
        TinySolver t1(1);
        TinySolver t2(2);
        TinySolver t3(3);
        assert(t1.solve_full() == 1);
        assert(t2.solve_full() == 6);
        assert(t3.solve_full() == 1008);

        const std::vector<u64> small_formula = compute_f_mod(6);
        for (int n = 1; n <= 6; ++n) {
            TinySolver t(n);
            const u64 brute_mod = mod_from_cpp_int(t.solve_full());
            assert(brute_mod == small_formula[static_cast<std::size_t>(n)]);
        }
    }

    const std::vector<u64> f = compute_f_mod(10'000);
    u64 answer = 0ULL;
    for (int n = 1; n <= 10'000; ++n) {
        answer += f[static_cast<std::size_t>(n)];
        if (answer >= kMod) {
            answer %= kMod;
        }
    }

    std::cout << answer << '\n';
    return 0;
}
