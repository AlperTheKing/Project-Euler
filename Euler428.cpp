#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

using int64 = long long;
using i128 = __int128_t;

namespace {

int64 isqrt_ll(int64 n) {
    int64 r = static_cast<int64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

class NecklaceCounter {
public:
    explicit NecklaceCounter(int64 max_n) {
        init_mu(max_n);
        reserve_caches();
    }

    int64 T(int64 n) {
        // k=3: (a-3b)(c-3b)=12 b^2, k=4: (a-b)(c-b)=2 b^2, k=6: (3a-b)(3c-b)=4 b^2.
        int64 s2 = S_D(n, 2);
        int64 s12 = S_D(n, 12);
        int64 s4 = S_D(n, 4);
        int64 s36_n3 = S_D(n / 3, 36);

        i128 c_not3 = (static_cast<i128>(s4) - s36_n3 - S_h(n)) / 2;

        i128 c_div3 = 0;
        int64 m = n / 3;
        int t = 1;
        while (m > 0) {
            i128 term = static_cast<i128>(S_D(m, 4)) - S_D(m / 3, 36);
            c_div3 += static_cast<i128>(2 * t - 1) * term;
            ++t;
            m /= 3;
        }

        i128 total = static_cast<i128>(s2) + s12 + c_not3 + c_div3;
        return static_cast<int64>(total);
    }

private:
    std::vector<int> mu;

    std::unordered_map<int64, int64> cache_c6;
    std::unordered_map<int64, int64> cache_a0;
    std::unordered_map<int64, int64> cache_chisq;
    std::unordered_map<int64, int64> cache_fg;
    std::unordered_map<int64, int64> cache_sh;
    std::unordered_map<int64, int64> cache_sd2;
    std::unordered_map<int64, int64> cache_sd4;
    std::unordered_map<int64, int64> cache_sd12;
    std::unordered_map<int64, int64> cache_sd36;

    void reserve_caches() {
        const size_t reserve = 1 << 18;
        cache_c6.reserve(reserve);
        cache_a0.reserve(reserve);
        cache_chisq.reserve(reserve);
        cache_fg.reserve(reserve);
        cache_sh.reserve(reserve);
        cache_sd2.reserve(reserve);
        cache_sd4.reserve(reserve);
        cache_sd12.reserve(reserve);
        cache_sd36.reserve(reserve);
    }

    void init_mu(int64 max_n) {
        int limit = static_cast<int>(isqrt_ll(max_n)) + 5;
        mu.assign(limit + 1, 0);
        std::vector<int> primes;
        std::vector<int> is_comp(limit + 1, 0);
        mu[1] = 1;
        for (int i = 2; i <= limit; ++i) {
            if (!is_comp[i]) {
                primes.push_back(i);
                mu[i] = -1;
            }
            for (int p : primes) {
                int64 v = static_cast<int64>(i) * p;
                if (v > limit) break;
                is_comp[static_cast<int>(v)] = 1;
                if (i % p == 0) {
                    mu[static_cast<int>(v)] = 0;
                    break;
                } else {
                    mu[static_cast<int>(v)] = -mu[i];
                }
            }
        }
    }

    int64 squarefree_count(int64 n) {
        if (n <= 0) return 0;
        int64 r = isqrt_ll(n);
        i128 sum = 0;
        for (int64 k = 1; k <= r; ++k) {
            sum += static_cast<i128>(mu[static_cast<size_t>(k)]) * (n / (k * k));
        }
        return static_cast<int64>(sum);
    }

    int64 S_c6(int64 n) {
        if (n <= 0) return 0;
        auto it = cache_c6.find(n);
        if (it != cache_c6.end()) return it->second;
        int64 res = squarefree_count(n) - S_c6(n / 2) - S_c6(n / 3) - S_c6(n / 6);
        cache_c6[n] = res;
        return res;
    }

    int64 A0(int64 n) {
        if (n <= 0) return 0;
        auto it = cache_a0.find(n);
        if (it != cache_a0.end()) return it->second;
        // Sum floor(n/d) over squarefree d coprime to 6.
        i128 res = 0;
        int64 l = 1;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            res += static_cast<i128>(v) * (S_c6(r) - S_c6(l - 1));
            l = r + 1;
        }
        int64 ans = static_cast<int64>(res);
        cache_a0[n] = ans;
        return ans;
    }

    int64 F_D(int64 n, int D) {
        if (n <= 0) return 0;
        if (D == 1) {
            return A0(n) + A0(n / 2) + A0(n / 3) + A0(n / 6);
        }
        if (D == 2) {
            return 2 * A0(n) + 2 * A0(n / 3);
        }
        if (D == 4) {
            return 3 * A0(n) + 3 * A0(n / 3) - A0(n / 2) - A0(n / 6);
        }
        if (D == 12) {
            return 6 * A0(n) - 2 * A0(n / 2);
        }
        if (D == 36) {
            return 9 * A0(n) - 3 * A0(n / 2) - 3 * A0(n / 3) + A0(n / 6);
        }
        return 0;
    }

    int64 S_D(int64 n, int D) {
        if (n <= 0) return 0;
        std::unordered_map<int64, int64>* cache = nullptr;
        if (D == 2) cache = &cache_sd2;
        else if (D == 4) cache = &cache_sd4;
        else if (D == 12) cache = &cache_sd12;
        else if (D == 36) cache = &cache_sd36;
        else return 0;

        auto it = cache->find(n);
        if (it != cache->end()) return it->second;

        i128 res = 0;
        int64 l = 1;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            res += static_cast<i128>(v) * (F_D(r, D) - F_D(l - 1, D));
            l = r + 1;
        }
        int64 ans = static_cast<int64>(res);
        (*cache)[n] = ans;
        return ans;
    }

    int64 sum_chi_prefix(int64 n) {
        if (n <= 0) return 0;
        return (n % 3 == 1) ? 1 : 0;
    }

    int64 S_chisq(int64 n) {
        if (n <= 0) return 0;
        auto it = cache_chisq.find(n);
        if (it != cache_chisq.end()) return it->second;
        // Sum_{m<=n} mu^2(m) * chi(m) via square divisors.
        int64 r = isqrt_ll(n);
        i128 sum = 0;
        for (int64 k = 1; k <= r; ++k) {
            if (k % 3 == 0) continue;
            sum += static_cast<i128>(mu[static_cast<size_t>(k)]) *
                   sum_chi_prefix(n / (k * k));
        }
        int64 ans = static_cast<int64>(sum);
        cache_chisq[n] = ans;
        return ans;
    }

    int64 F_G(int64 n) {
        if (n <= 0) return 0;
        auto it = cache_fg.find(n);
        if (it != cache_fg.end()) return it->second;
        i128 res = 0;
        int64 l = 1;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            res += static_cast<i128>(v) * (S_chisq(r) - S_chisq(l - 1));
            l = r + 1;
        }
        int64 ans = static_cast<int64>(res);
        cache_fg[n] = ans;
        return ans;
    }

    int64 F1(int64 n) {
        if (n <= 0) return 0;
        return F_G(n) - F_G(n / 3);
    }

    int64 S_h(int64 n) {
        if (n <= 0) return 0;
        auto it = cache_sh.find(n);
        if (it != cache_sh.end()) return it->second;
        // Convolution with the mod-3 character yields a simple filter on floor(n/d) mod 3.
        i128 res = 0;
        int64 l = 1;
        while (l <= n) {
            int64 v = n / l;
            int64 r = n / v;
            if (v % 3 == 1) {
                res += static_cast<i128>(F1(r) - F1(l - 1));
            }
            l = r + 1;
        }
        int64 ans = static_cast<int64>(res);
        cache_sh[n] = ans;
        return ans;
    }
};

}  // namespace

int main() {
    const int64 target = 1000000000LL;
    NecklaceCounter solver(target);

    struct Check { int64 n; int64 expected; };
    const Check checks[] = {
        {1, 9},
        {20, 732},
        {3000, 438106},
    };

    for (const auto& chk : checks) {
        int64 got = solver.T(chk.n);
        if (got != chk.expected) {
            std::cerr << "Validation failed for n=" << chk.n
                      << ": got " << got << ", expected " << chk.expected << '\n';
            return 1;
        }
        std::cout << "Validation passed for n=" << chk.n << '\n';
    }

    std::cout << solver.T(target) << '\n';
    return 0;
}
