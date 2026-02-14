#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int MOD = 1'000'000'007;
constexpr u64 MAIN_N = 1'000'000'000'000ULL;
constexpr int K_MAX = 50;

int mod_pow(long long a, long long e) {
    long long r = 1 % MOD;
    a %= MOD;
    while (e > 0) {
        if (e & 1LL) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1LL;
    }
    return static_cast<int>(r);
}

int add_mod(int a, int b) {
    int s = a + b;
    if (s >= MOD) s -= MOD;
    return s;
}

int sub_mod(int a, int b) {
    int s = a - b;
    if (s < 0) s += MOD;
    return s;
}

int mul_mod(long long a, long long b) {
    return static_cast<int>((static_cast<u128>(a) * b) % MOD);
}

std::vector<int> primes_up_to(int n) {
    if (n < 2) return {};
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(n + 1), 1U);
    is_prime[0] = is_prime[1] = 0U;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) continue;
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = 0U;
        }
    }
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

struct FaulhaberData {
    int K = 0;
    std::vector<std::vector<int>> coeff;

    explicit FaulhaberData(int max_k) : K(max_k), coeff(static_cast<std::size_t>(max_k + 1)) {
        std::vector<int> fact(static_cast<std::size_t>(K + 2), 1);
        std::vector<int> inv_fact(static_cast<std::size_t>(K + 2), 1);
        std::vector<int> inv(static_cast<std::size_t>(K + 3), 0);

        for (int i = 1; i <= K + 1; ++i) {
            fact[static_cast<std::size_t>(i)] = mul_mod(fact[static_cast<std::size_t>(i - 1)], i);
        }
        inv_fact[static_cast<std::size_t>(K + 1)] = mod_pow(fact[static_cast<std::size_t>(K + 1)], MOD - 2);
        for (int i = K + 1; i >= 1; --i) {
            inv_fact[static_cast<std::size_t>(i - 1)] =
                mul_mod(inv_fact[static_cast<std::size_t>(i)], i);
        }
        for (int i = 1; i <= K + 2; ++i) {
            inv[static_cast<std::size_t>(i)] = mod_pow(i, MOD - 2);
        }

        auto ncr = [&](int n, int r) -> int {
            if (r < 0 || r > n) return 0;
            return mul_mod(mul_mod(fact[static_cast<std::size_t>(n)], inv_fact[static_cast<std::size_t>(r)]),
                           inv_fact[static_cast<std::size_t>(n - r)]);
        };

        std::vector<int> B(static_cast<std::size_t>(K + 1), 0);
        B[0] = 1;
        if (K >= 1) B[1] = (MOD + 1) / 2;

        for (int k = 2; k <= K; k += 2) {
            int s = 0;
            for (int i = 0; i < k; ++i) {
                int t = mul_mod(ncr(k, i), B[static_cast<std::size_t>(i)]);
                t = mul_mod(t, inv[static_cast<std::size_t>(k - i + 1)]);
                s = add_mod(s, t);
            }
            B[static_cast<std::size_t>(k)] = sub_mod(1, s);
        }

        for (int k = 1; k <= K; ++k) {
            std::vector<int> v(static_cast<std::size_t>(k + 2), 0);
            for (int i = 1; i < k; ++i) {
                int t = B[static_cast<std::size_t>(k + 1 - i)];
                t = mul_mod(t, ncr(k, i));
                t = mul_mod(t, inv[static_cast<std::size_t>(k + 1 - i)]);
                v[static_cast<std::size_t>(i)] = t;
            }
            v[static_cast<std::size_t>(k)] = (MOD + 1) / 2;
            v[static_cast<std::size_t>(k + 1)] = inv[static_cast<std::size_t>(k + 1)];
            coeff[static_cast<std::size_t>(k)] = std::move(v);
        }
    }

    int eval(u64 n, int k) const {
        const auto& v = coeff[static_cast<std::size_t>(k)];
        int n_mod = static_cast<int>(n % static_cast<u64>(MOD));
        int m = 1;
        int acc = 0;
        for (int c : v) {
            acc = add_mod(acc, mul_mod(m, c));
            m = mul_mod(m, n_mod);
        }
        return acc;
    }
};

struct ValData {
    std::vector<int> primes;
    std::vector<std::vector<u64>> V;
};

void collect_powerful(std::vector<u64>& out, u64 m, int i, u64 n, const std::vector<int>& primes) {
    out.push_back(m);

    if (static_cast<u128>(primes[static_cast<std::size_t>(i)]) * m > n) {
        return;
    }

    collect_powerful(out, m * static_cast<u64>(primes[static_cast<std::size_t>(i)]), i, n, primes);

    const int lP = static_cast<int>(primes.size());
    for (int j = i + 1; j < lP; ++j) {
        const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(j)]);
        const u128 mm = static_cast<u128>(m) * p * p;
        if (mm > n) {
            return;
        }
        collect_powerful(out, static_cast<u64>(mm), j, n, primes);
    }
}

ValData build_val_data(u64 n) {
    const int lim = static_cast<int>(std::sqrt(static_cast<long double>(n)));
    std::vector<int> P = primes_up_to(lim);
    const int lP = static_cast<int>(P.size());

    std::vector<std::vector<u64>> V(static_cast<std::size_t>(lP + 1));
    V[static_cast<std::size_t>(lP)] = {n};

    std::unordered_set<u64> S;
    S.reserve(1 << 16);

    for (int i = lP - 1; i >= 0; --i) {
        std::vector<u64> arr;
        arr.reserve(256);
        arr.push_back(1ULL);

        const u64 p = static_cast<u64>(P[static_cast<std::size_t>(i)]);
        collect_powerful(arr, p * p, i, n, P);

        for (u64 x : arr) {
            S.insert(n / x);
        }

        std::vector<u64> cur;
        cur.reserve(S.size());
        for (u64 x : S) {
            cur.push_back(x);
        }
        std::sort(cur.begin(), cur.end(), std::greater<u64>());
        V[static_cast<std::size_t>(i)] = std::move(cur);
    }

    return {std::move(P), std::move(V)};
}

int solve_total(u64 n, int K) {
    const FaulhaberData F(K);
    const ValData data = build_val_data(n);

    const auto& P = data.primes;
    const auto& V = data.V;
    const int lP = static_cast<int>(P.size());

    const auto& V0 = V[0];
    std::unordered_map<u64, int> idx;
    idx.reserve(V0.size() * 2 + 16);
    for (int i = 0; i < static_cast<int>(V0.size()); ++i) {
        idx.emplace(V0[static_cast<std::size_t>(i)], i);
    }

    std::vector<std::vector<int>> Vidx(static_cast<std::size_t>(lP + 1));
    for (int i = 0; i <= lP; ++i) {
        const auto& vec = V[static_cast<std::size_t>(i)];
        auto& idv = Vidx[static_cast<std::size_t>(i)];
        idv.reserve(vec.size());
        for (u64 x : vec) {
            idv.push_back(idx[x]);
        }
    }

    int total = 0;
    for (int k = 1; k <= K; ++k) {
        std::vector<int> Svals(V0.size(), 0);
        for (int id : Vidx[0]) {
            Svals[static_cast<std::size_t>(id)] = F.eval(V0[static_cast<std::size_t>(id)], k);
        }

        for (int i = 0; i < lP; ++i) {
            const u64 p = static_cast<u64>(P[static_cast<std::size_t>(i)]);
            const u64 p2 = p * p;
            const int pk = mod_pow(static_cast<int>(p % MOD), k);
            const int t = mul_mod(pk, sub_mod(pk, 1));

            const auto& vec = V[static_cast<std::size_t>(i + 1)];
            const auto& idv = Vidx[static_cast<std::size_t>(i + 1)];

            for (std::size_t pos = 0; pos < vec.size(); ++pos) {
                const u64 x = vec[pos];
                if (x < p2) break;

                const int idx_x = idv[pos];
                int cur = Svals[static_cast<std::size_t>(idx_x)];

                u64 pp = p2;
                while (pp <= x) {
                    const auto it = idx.find(x / pp);
                    cur = sub_mod(cur, mul_mod(t, Svals[static_cast<std::size_t>(it->second)]));
                    if (pp > x / p) break;
                    pp *= p;
                }
                Svals[static_cast<std::size_t>(idx_x)] = cur;
            }
        }

        total = add_mod(total, Svals[static_cast<std::size_t>(idx[n])]);
    }

    return total;
}

bool run_validations() {
    const int s1_10 = solve_total(10ULL, 1);
    if (s1_10 != 41) {
        std::cerr << "Validation failed: S_1(10)=" << s1_10 << "\n";
        return false;
    }

    const int s1_100 = solve_total(100ULL, 1);
    if (s1_100 != 3512) {
        std::cerr << "Validation failed: S_1(100)=" << s1_100 << "\n";
        return false;
    }

    const int sum2_100 = solve_total(100ULL, 2);
    const int s2_100 = sub_mod(sum2_100, s1_100);
    if (s2_100 != 208090) {
        std::cerr << "Validation failed: S_2(100)=" << s2_100 << "\n";
        return false;
    }

    const int s1_10000 = solve_total(10000ULL, 1);
    if (s1_10000 != 35252550) {
        std::cerr << "Validation failed: S_1(10000)=" << s1_10000 << "\n";
        return false;
    }

    const int sum3_1e8 = solve_total(100000000ULL, 3);
    if (sum3_1e8 != 338787512) {
        std::cerr << "Validation failed: sum_{k<=3} S_k(1e8)=" << sum3_1e8 << "\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool validate = true;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--no-validate") {
            validate = false;
        }
    }

    if (validate && !run_validations()) return 1;
    std::cout << solve_total(MAIN_N, K_MAX) << '\n';
    return 0;
}
