#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using i64 = std::int64_t;

static constexpr int MOD = 1'000'000'007;

static int mod_pow(i64 a, i64 e) {
    i64 r = 1;
    a %= MOD;
    while (e > 0) {
        if (e & 1LL) {
            r = (r * a) % MOD;
        }
        a = (a * a) % MOD;
        e >>= 1LL;
    }
    return static_cast<int>(r);
}

static inline int mod_mul(int a, int b) {
    return static_cast<int>((static_cast<i64>(a) * static_cast<i64>(b)) % MOD);
}

static int solve_fast(int N) {
    std::vector<int> lp(N + 1, 0);
    std::vector<int> primes;
    primes.reserve(N / 10);

    std::vector<std::int8_t> mu(N + 1, 0);
    mu[1] = 1;

    for (int i = 2; i <= N; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            const i64 v = static_cast<i64>(p) * static_cast<i64>(i);
            if (v > N) {
                break;
            }
            lp[static_cast<int>(v)] = p;
            if (p == lp[i]) {
                mu[static_cast<int>(v)] = 0;
                break;
            }
            mu[static_cast<int>(v)] = static_cast<std::int8_t>(-mu[i]);
        }
    }

    lp.clear();
    lp.shrink_to_fit();

    std::vector<int> mertens_mod(N + 1, 0);
    for (int i = 1; i <= N; ++i) {
        int v = mertens_mod[i - 1] + static_cast<int>(mu[i]);
        if (v >= MOD) {
            v -= MOD;
        }
        if (v < 0) {
            v += MOD;
        }
        mertens_mod[i] = v;
    }

    std::vector<int> base(N + 1, 0);
    int p2 = 1;
    for (int i = 1; i <= N; ++i) {
        p2 = mod_mul(p2, 2);
        int v = p2 - 1;
        if (v < 0) {
            v += MOD;
        }
        base[i] = v;
    }

    std::vector<int> inv_base(N + 1, 0);
    inv_base[0] = 1;
    for (int i = 1; i <= N; ++i) {
        inv_base[i] = mod_mul(inv_base[i - 1], base[i]);
    }

    int inv_pref = mod_pow(inv_base[N], MOD - 2);
    for (int i = N; i >= 1; --i) {
        inv_base[i] = mod_mul(inv_pref, inv_base[i - 1]);
        inv_pref = mod_mul(inv_pref, base[i]);
    }

    std::vector<int> phi2(N + 1, 1);
    for (int k = 1; k <= N; ++k) {
        const int muk = static_cast<int>(mu[k]);
        if (muk == 0) {
            continue;
        }

        const int limit = N / k;
        int n = k;
        if (muk > 0) {
            for (int m = 1; m <= limit; ++m, n += k) {
                phi2[n] = mod_mul(phi2[n], base[m]);
            }
        } else {
            for (int m = 1; m <= limit; ++m, n += k) {
                phi2[n] = mod_mul(phi2[n], inv_base[m]);
            }
        }
    }

    base.clear();
    base.shrink_to_fit();
    inv_base.clear();
    inv_base.shrink_to_fit();

    std::vector<int> F(N + 1, 1);
    for (int d = 1; d <= N; ++d) {
        int factor = phi2[d] + 1;
        if (factor >= MOD) {
            factor -= MOD;
        }
        for (int n = d; n <= N; n += d) {
            F[n] = mod_mul(F[n], factor);
        }
    }

    int ans = 0;
    for (int n = 1; n <= N; ++n) {
        ans = (ans + static_cast<i64>(F[n]) * static_cast<i64>(mertens_mod[N / n])) % MOD;
    }

    return ans;
}

static i64 brute_q_small(int N) {
    std::vector<i64> cyclo(N + 1, 0);
    cyclo[1] = 1;

    for (int n = 2; n <= N; ++n) {
        i64 v = (1LL << n) - 1;
        for (int d = 1; d < n; ++d) {
            if (n % d == 0) {
                v /= cyclo[d];
            }
        }
        cyclo[n] = v;
    }

    i64 q = 0;
    for (int n = 1; n <= N; ++n) {
        std::vector<int> divs;
        for (int d = 1; d <= n; ++d) {
            if (n % d == 0) {
                divs.push_back(d);
            }
        }

        i64 pn = 0;
        const int z = static_cast<int>(divs.size());
        for (int mask = 0; mask < (1 << z); ++mask) {
            int l = 1;
            i64 prod = 1;
            for (int i = 0; i < z; ++i) {
                if ((mask >> i) & 1) {
                    l = std::lcm(l, divs[i]);
                    prod *= cyclo[divs[i]];
                }
            }
            if (l == n) {
                pn += prod;
            }
        }
        q += pn;
    }

    return q;
}

int main() {
    assert(solve_fast(10) == 5'598);
    assert(solve_fast(8) == brute_q_small(8));

    std::cout << solve_fast(10'000'000) << '\n';
    return 0;
}
