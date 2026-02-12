#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u16 = std::uint16_t;
using u128 = unsigned __int128;

static constexpr i64 MOD = 1'000'000'007LL;

static i64 mod_norm(i64 x) {
    x %= MOD;
    if (x < 0) x += MOD;
    return x;
}

static i64 mod_mul(i64 a, i64 b) {
    return static_cast<i64>((static_cast<u128>(mod_norm(a)) * static_cast<u128>(mod_norm(b))) % MOD);
}

static std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int p = 2; p * p <= n; ++p) {
        if (!is_prime[p]) continue;
        for (int q = p * p; q <= n; q += p) is_prime[q] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) if (is_prime[i]) primes.push_back(i);
    return primes;
}

struct Engine {
    int N;
    std::vector<int> small_primes;
    std::vector<int> large_primes;
    std::vector<int> max_exp;
    std::vector<int> radix;
    std::vector<int> stride;
    int states = 1;
    std::vector<std::vector<i64>> p_pow;
    std::vector<std::pair<int, int>> small_numbers;

    explicit Engine(int n) : N(n) {
        auto primes = primes_up_to(N);
        int root = 1;
        while ((root + 1) * (root + 1) <= N) ++root;

        for (int p : primes) {
            if (p <= root) small_primes.push_back(p);
            else large_primes.push_back(p);
        }

        int k = static_cast<int>(small_primes.size());
        max_exp.assign(k, 0);
        radix.assign(k, 0);

        for (int i = 0; i < k; ++i) {
            int p = small_primes[i];
            int e = 0;
            i64 v = 1;
            while (v * p <= N) {
                v *= p;
                ++e;
            }
            max_exp[i] = e;
            radix[i] = e + 1;
            states *= radix[i];
        }

        stride.assign(k, 1);
        for (int i = k - 2; i >= 0; --i) stride[i] = stride[i + 1] * radix[i + 1];

        p_pow.resize(k);
        for (int i = 0; i < k; ++i) {
            int emax = max_exp[i];
            p_pow[i].assign(emax + 1, 1);
            for (int e = 1; e <= emax; ++e) p_pow[i][e] = mod_mul(p_pow[i][e - 1], small_primes[i]);
        }

        std::vector<int> exps(k, 0);
        gen_small_numbers(0, 1, 0, exps);
        std::sort(small_numbers.begin(), small_numbers.end());
    }

    void gen_small_numbers(int i, i64 value, int idx, std::vector<int>& exps) {
        if (i == static_cast<int>(small_primes.size())) {
            small_numbers.emplace_back(static_cast<int>(value), idx);
            return;
        }

        i64 v = value;
        int p = small_primes[i];
        int base_idx = idx * radix[i];
        for (int e = 0; e <= max_exp[i]; ++e) {
            exps[i] = e;
            gen_small_numbers(i + 1, v, base_idx + e, exps);
            if (e == max_exp[i]) break;
            if (v > N / p) break;
            v *= p;
        }
    }

    void prefix_transform(std::vector<int>& arr) const {
        int k = static_cast<int>(small_primes.size());
        for (int d = 0; d < k; ++d) {
            int st = stride[d];
            int block = radix[d] * st;
            for (int base = 0; base < states; base += block) {
                for (int off = 0; off < st; ++off) {
                    int idx = base + off;
                    for (int e = 1; e < radix[d]; ++e) {
                        idx += st;
                        arr[idx] += arr[idx - st];
                    }
                }
            }
        }
    }

    std::vector<u16> counts_for_limit(int limit) const {
        std::vector<int> arr(states, 0);
        for (const auto& [val, idx] : small_numbers) {
            if (val > limit) break;
            arr[idx] = 1;
        }
        prefix_transform(arr);

        std::vector<u16> out(states);
        for (int i = 0; i < states; ++i) out[i] = static_cast<u16>(arr[i]);
        return out;
    }

    static int floor_div(int a, int b) {
        return a / b;
    }

    i64 solve_mod() {
        std::vector<int> needed_limits;
        needed_limits.push_back(N);

        std::unordered_map<int, std::vector<int>> groups;
        groups.reserve(64);
        for (int p : large_primes) {
            int L = floor_div(N, p);
            groups[L].push_back(p);
        }
        for (auto& kv : groups) needed_limits.push_back(kv.first);

        std::sort(needed_limits.begin(), needed_limits.end());
        needed_limits.erase(std::unique(needed_limits.begin(), needed_limits.end()), needed_limits.end());

        std::unordered_map<int, std::vector<u16>> counts;
        counts.reserve(needed_limits.size() * 2);
        for (int L : needed_limits) counts[L] = counts_for_limit(L);

        std::vector<i64> pow2(N + 1, 1);
        for (int i = 1; i <= N; ++i) pow2[i] = (pow2[i - 1] * 2) % MOD;

        std::unordered_map<int, std::vector<i64>> group_value;
        group_value.reserve(groups.size() * 2);
        for (const auto& kv : groups) {
            int L = kv.first;
            const auto& primes = kv.second;
            std::vector<i64> table(N + 1, 1);
            for (int t = 0; t <= N; ++t) {
                i64 v = 1;
                for (int p : primes) {
                    i64 term = mod_norm(1 - p + mod_mul(p, pow2[t]));
                    v = mod_mul(v, term);
                }
                table[t] = v;
            }
            group_value[L] = std::move(table);
        }

        i64 ans = 0;
        int k = static_cast<int>(small_primes.size());

        for (int idx = 0; idx < states; ++idx) {
            int x = idx;
            i64 d = 1;
            i64 cut = 1;

            for (int i = k - 1; i >= 0; --i) {
                int e = x % radix[i];
                x /= radix[i];

                d = mod_mul(d, p_pow[i][e]);
                if (e < max_exp[i]) cut = mod_mul(cut, mod_norm(1 - small_primes[i]));
            }

            i64 term = mod_mul(d, cut);
            int t0 = counts[N][idx];
            term = mod_mul(term, pow2[t0]);

            for (const auto& kv : group_value) {
                int L = kv.first;
                int t = counts[L][idx];
                term = mod_mul(term, kv.second[t]);
            }

            ans += term;
            if (ans >= MOD) ans -= MOD;
        }

        return ans;
    }
};

static i64 G_mod(int N) {
    Engine e(N);
    return e.solve_mod();
}

int main() {
    assert(G_mod(5) == 528);
    assert(G_mod(20) == 108589719);
    std::cout << G_mod(800) << '\n';
    return 0;
}
