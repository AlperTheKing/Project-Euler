#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = long long;
using boost::multiprecision::cpp_int;

static u32 mul_mod(u32 a, u32 b, u32 mod) {
    return static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(b)) % static_cast<u64>(mod));
}

static u32 mod_pow(u32 a, u64 e, u32 mod) {
    u32 r = 1;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod(r, a, mod);
        }
        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static std::vector<int> primes_upto(int n) {
    std::vector<bool> sieve(n + 1, true);
    sieve[0] = sieve[1] = false;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!sieve[i]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            sieve[j] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (sieve[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

static void enumerate_min_costs(
    const std::vector<int>& primes,
    int idx,
    int cost,
    u32 residue,
    int limit,
    u32 mod,
    std::unordered_map<u32, int>& best
) {
    if (idx == static_cast<int>(primes.size())) {
        auto it = best.find(residue);
        if (it == best.end() || cost < it->second) {
            best[residue] = cost;
        }
        return;
    }

    const int q = primes[idx];
    const int w = q - 1;
    const int max_e = (limit - cost) / w;

    u32 r = residue;
    for (int e = 0; e <= max_e; ++e) {
        enumerate_min_costs(primes, idx + 1, cost + e * w, r, limit, mod, best);
        r = mul_mod(r, static_cast<u32>(q), mod);
    }
}

static bool find_combo_exact(
    const std::vector<int>& primes,
    int idx,
    int rem_cost,
    u32 residue,
    u32 target,
    u32 mod,
    std::vector<int>& exps
) {
    if (idx == static_cast<int>(primes.size())) {
        return rem_cost == 0 && residue == target;
    }

    const int q = primes[idx];
    const int w = q - 1;
    const int max_e = rem_cost / w;

    u32 r = residue;
    for (int e = 0; e <= max_e; ++e) {
        exps[idx] = e;
        if (find_combo_exact(primes, idx + 1, rem_cost - e * w, r, target, mod, exps)) {
            return true;
        }
        r = mul_mod(r, static_cast<u32>(q), mod);
    }

    return false;
}

struct SolveResult {
    int extra_cost;
    cpp_int product;
};

static SolveResult solve_prime(u32 p) {
    int limit = 80;

    while (true) {
        const int pmax = std::min<int>(static_cast<int>(p) - 1, limit + 1);
        auto primes = primes_upto(pmax);
        const int split = std::min<int>(4, static_cast<int>(primes.size()));

        std::vector<int> pa(primes.begin(), primes.begin() + split);
        std::vector<int> pb(primes.begin() + split, primes.end());

        std::unordered_map<u32, int> bestA;
        std::unordered_map<u32, int> bestB;

        if (p > 1'000'000U) {
            bestA.reserve(5'000'000);
            bestB.reserve(2'000'000);
        }

        enumerate_min_costs(pa, 0, 0, 1, limit, p, bestA);
        enumerate_min_costs(pb, 0, 0, 1, limit, p, bestB);

        const u32 target = p - 1;

        std::vector<u32> keys;
        std::vector<int> costsA;
        keys.reserve(bestA.size());
        costsA.reserve(bestA.size());
        for (const auto& [res, c] : bestA) {
            keys.push_back(res);
            costsA.push_back(c);
        }

        std::vector<u32> prefix(keys.size() + 1, 1);
        for (std::size_t i = 0; i < keys.size(); ++i) {
            prefix[i + 1] = mul_mod(prefix[i], keys[i], p);
        }
        u32 inv_total = mod_pow(prefix.back(), static_cast<u64>(p) - 2ULL, p);

        std::vector<u32> invs(keys.size(), 1);
        for (std::size_t i = keys.size(); i-- > 0;) {
            invs[i] = mul_mod(inv_total, prefix[i], p);
            inv_total = mul_mod(inv_total, keys[i], p);
        }

        int best_total = std::numeric_limits<int>::max();
        int best_ca = -1;
        int best_cb = -1;
        u32 best_ra = 0;
        u32 best_rb = 0;

        for (std::size_t i = 0; i < keys.size(); ++i) {
            const u32 ra = keys[i];
            const int ca = costsA[i];
            const u32 need = mul_mod(target, invs[i], p);
            const auto it = bestB.find(need);
            if (it == bestB.end()) {
                continue;
            }

            const int cb = it->second;
            const int total = ca + cb;
            if (total < best_total) {
                best_total = total;
                best_ca = ca;
                best_cb = cb;
                best_ra = ra;
                best_rb = need;
            }
        }

        if (best_total != std::numeric_limits<int>::max() && best_total <= limit) {
            std::vector<int> expsA(pa.size(), 0);
            std::vector<int> expsB(pb.size(), 0);

            const bool okA = find_combo_exact(pa, 0, best_ca, 1, best_ra, p, expsA);
            const bool okB = find_combo_exact(pb, 0, best_cb, 1, best_rb, p, expsB);
            if (!okA || !okB) {
                throw std::runtime_error("reconstruction failed");
            }

            cpp_int prod = 1;
            int cost_check = 0;
            u32 residue_check = 1;

            for (int i = 0; i < static_cast<int>(pa.size()); ++i) {
                const int q = pa[i];
                const int e = expsA[i];
                for (int j = 0; j < e; ++j) {
                    prod *= q;
                    residue_check = mul_mod(residue_check, static_cast<u32>(q), p);
                    cost_check += q - 1;
                }
            }
            for (int i = 0; i < static_cast<int>(pb.size()); ++i) {
                const int q = pb[i];
                const int e = expsB[i];
                for (int j = 0; j < e; ++j) {
                    prod *= q;
                    residue_check = mul_mod(residue_check, static_cast<u32>(q), p);
                    cost_check += q - 1;
                }
            }

            if (cost_check != best_total || residue_check != target) {
                throw std::runtime_error("consistency check failed");
            }

            return {best_total, prod};
        }

        limit += 40;
        if (limit > 600) {
            throw std::runtime_error("failed to find solution within limit");
        }
    }
}

static int brute_min_total_cost(int p) {
    const int n = p - 1;
    const int full = (1 << n) - 1;

    std::vector<int> memo(1 << n, -1);

    std::function<int(int)> dfs = [&](int mask) -> int {
        if (mask == 0) {
            return 0;
        }
        int& ret = memo[mask];
        if (ret != -1) {
            return ret;
        }

        int i = 0;
        while (((mask >> i) & 1) == 0) {
            ++i;
        }

        int best = std::numeric_limits<int>::max();
        const int rem = mask ^ (1 << i);
        for (int j = i + 1; j < n; ++j) {
            if (((rem >> j) & 1) == 0) {
                continue;
            }
            const int c = static_cast<int>((1LL * (i + 1) * (j + 1)) % p);
            best = std::min(best, c + dfs(rem ^ (1 << j)));
        }

        ret = best;
        return ret;
    };

    return dfs(full);
}

int main() {
    {
        const auto r5 = solve_prime(5);
        assert(r5.extra_cost == 2);
        assert(r5.product == 4);
    }

    for (int p : {7, 11, 13, 17, 19}) {
        const auto r = solve_prime(static_cast<u32>(p));
        const int brute = brute_min_total_cost(p);
        assert(brute == r.extra_cost + (p - 1) / 2);
    }

    const auto ans = solve_prime(2'000'000'011u);
    std::cout << ans.product << '\n';
    return 0;
}
