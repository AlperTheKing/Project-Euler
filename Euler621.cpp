#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;
using i128 = __int128_t;
using u128 = __uint128_t;

static u64 mod_mul(u64 a, u64 b, u64 mod) { return (u128)a * b % mod; }

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    while (e) {
        if (e & 1) r = mod_mul(r, a, mod);
        a = mod_mul(a, a, mod);
        e >>= 1;
    }
    return r;
}

static bool is_prime(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return n == p;
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    auto witness = [&](u64 a) -> bool {
        if (a % n == 0) return false;
        u64 x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) return false;
        for (int i = 1; i < s; ++i) {
            x = mod_mul(x, x, n);
            if (x == n - 1) return false;
        }
        return true;
    };

    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (witness(a)) return false;
    }
    return true;
}

static u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    static std::mt19937_64 rng(1234567);
    std::uniform_int_distribution<u64> dist(0, n - 1);

    while (true) {
        u64 c = dist(rng) % n;
        u64 x = dist(rng) % n;
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) { return (mod_mul(v, v, n) + c) % n; };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            u64 diff = (x > y) ? (x - y) : (y - x);
            d = std::gcd(diff, n);
        }
        if (d != n) return d;
    }
}

static void factor_rec(u64 n, std::map<u64, int> &out) {
    if (n == 1) return;
    if (is_prime(n)) {
        out[n]++;
        return;
    }
    u64 d = pollard_rho(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}

static u64 tonelli_shanks(u64 n, u64 p) {
    if (n == 0) return 0;
    if (p == 2) return n;

    if (mod_pow(n, (p - 1) / 2, p) != 1) return 0;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);

    u64 q = p - 1;
    int s = 0;
    while ((q & 1) == 0) {
        q >>= 1;
        ++s;
    }

    u64 z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) ++z;

    u64 c = mod_pow(z, q, p);
    u64 r = mod_pow(n, (q + 1) / 2, p);
    u64 t = mod_pow(n, q, p);
    int m = s;

    while (t != 1) {
        int i = 1;
        u64 tt = mod_mul(t, t, p);
        while (i < m && tt != 1) {
            tt = mod_mul(tt, tt, p);
            ++i;
        }
        u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        r = mod_mul(r, b, p);
        u64 b2 = mod_mul(b, b, p);
        t = mod_mul(t, b2, p);
        c = b2;
        m = i;
    }
    return r;
}

static u64 mod_norm_i64(i64 a, u64 mod) {
    i128 r = (i128)a % (i128)mod;
    if (r < 0) r += (i128)mod;
    return (u64)r;
}

static u64 inv_mod_u64(u64 a, u64 mod) {
    i128 t = 0, newt = 1;
    i128 r = (i128)mod, newr = (i128)a;
    while (newr != 0) {
        u64 q = (u64)(r / newr);
        i128 tmp = t - (i128)q * newt;
        t = newt;
        newt = tmp;
        tmp = r - (i128)q * newr;
        r = newr;
        newr = tmp;
    }
    assert(r == 1);
    if (t < 0) t += (i128)mod;
    return (u64)t;
}

static std::vector<u64> roots_mod_prime_power(i64 D, u64 p, int e) {
    u64 pe = 1;
    for (int i = 0; i < e; ++i) pe *= p;

    if (mod_norm_i64(D, p) == 0) {
        if (e == 1) return {0};
        return {};
    }

    u64 n = mod_norm_i64(D, p);
    u64 r = tonelli_shanks(n, p);
    if (r == 0) return {};

    u64 pk = p;
    for (int k = 1; k < e; ++k) {
        u64 pk1 = pk * p;
        u64 Dmod = mod_norm_i64(D, pk1);
        u64 r2 = (u128)r * r % pk1;
        u64 diff = (Dmod + pk1 - r2) % pk1;
        u64 delta = diff / pk;

        u64 inv = inv_mod_u64((2 * (r % p)) % p, p);
        u64 t = (u128)delta * inv % p;
        r += t * pk;
        pk = pk1;
    }

    r %= pe;
    u64 r2 = (pe - r) % pe;
    if (r2 == r) return {r};
    return {r, r2};
}

static std::vector<int> sieve_spf(int n) {
    std::vector<int> spf(n + 1, 0), primes;
    primes.reserve((size_t)n / 10);
    for (int i = 2; i <= n; ++i) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (int p : primes) {
            long long v = 1LL * p * i;
            if (v > n) break;
            spf[(int)v] = p;
            if (p == spf[i]) break;
        }
    }
    return spf;
}

static std::vector<std::pair<int, int>> factorize_small(int x, const std::vector<int> &spf) {
    std::vector<std::pair<int, int>> f;
    while (x > 1) {
        int p = spf[x];
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        f.push_back({p, e});
    }
    return f;
}

static i64 class_number_odd_fundamental(i64 D) {
    assert(D < 0);
    assert((D & 3LL) == 1); // D ≡ 1 (mod 4)
    const u64 absD = (u64)(-D);

    i64 maxA = (i64)std::sqrt((long double)absD / 3.0L);
    while (3 * (i128)(maxA + 1) * (maxA + 1) <= (i128)absD) ++maxA;
    while (3 * (i128)maxA * maxA > (i128)absD) --maxA;

    const int lim = (int)maxA;
    const std::vector<int> spf = sieve_spf(lim);

    i64 h = 0;
    for (int a = 1; a <= lim; a += 2) {
        const auto fac = factorize_small(a, spf);
        u64 mod = 1;
        std::vector<u64> roots = {0};
        bool ok = true;

        for (auto [p_i, e] : fac) {
            u64 p = (u64)p_i;
            u64 pe = 1;
            for (int i = 0; i < e; ++i) pe *= p;

            std::vector<u64> rpe = roots_mod_prime_power(D, p, e);
            if (rpe.empty()) {
                ok = false;
                break;
            }

            u64 inv = inv_mod_u64(mod % pe, pe);
            std::vector<u64> next;
            next.reserve(roots.size() * rpe.size());

            for (u64 x : roots) {
                u64 xpe = x % pe;
                for (u64 y : rpe) {
                    u64 t = (y + pe - xpe) % pe;
                    t = (u128)t * inv % pe;
                    next.push_back(x + mod * t);
                }
            }

            mod *= pe;
            roots.swap(next);
        }

        if (!ok) continue;
        assert(mod == (u64)a);

        for (u64 r : roots) {
            u64 bmod = r;
            if ((bmod & 1ULL) == 0) bmod += (u64)a;
            bmod %= 2ULL * (u64)a;

            i64 b = (i64)bmod;
            if (b > a) b -= 2LL * a;
            if (b == -a) b = a;

            const i128 num = (i128)b * b - (i128)D;
            const i128 den = 4 * (i128)a;
            if (num % den != 0) continue;
            const i64 c = (i64)(num / den);
            if (a > c) continue;
            if ((a == c || std::llabs(b) == a) && b < 0) continue;

            ++h;
        }
    }
    return h;
}

static int jacobi_symbol(i64 a, i64 n) {
    assert(n > 0 && (n & 1LL) == 1);
    a %= n;
    if (a < 0) a += n;
    int s = 1;
    while (a != 0) {
        while ((a & 1LL) == 0) {
            a >>= 1;
            i64 r = n & 7LL;
            if (r == 3 || r == 5) s = -s;
        }
        std::swap(a, n);
        if ((a & 3LL) == 3 && (n & 3LL) == 3) s = -s;
        a %= n;
    }
    return (n == 1) ? s : 0;
}

static int kronecker_D_over_2(i64 D) {
    assert(D & 1LL);
    i64 r = D & 7LL;
    if (r == 1 || r == 7) return 1;
    if (r == 3 || r == 5) return -1;
    assert(false);
    return 0;
}

static u64 sigma_prime_power(u64 p, int e) {
    u128 pe1 = 1;
    for (int i = 0; i < e + 1; ++i) pe1 *= p;
    return (u64)((pe1 - 1) / (p - 1));
}

static u64 r3_sum_of_three_squares(u64 N) {
    std::map<u64, int> pf;
    factor_rec(N, pf);

    u64 m = 1;
    std::map<u64, int> f_pf;
    for (auto [p, e] : pf) {
        if (e & 1) m *= p;
        int k = e / 2;
        if (k > 0) f_pf[p] = k;
    }

    const i64 D = -(i64)m;
    const i64 h = class_number_odd_fundamental(D);
    const int w = (D == -3) ? 6 : ((D == -4) ? 4 : 2);
    const int one_minus = 1 - kronecker_D_over_2(D);

    std::vector<u64> primes;
    primes.reserve(f_pf.size());
    std::vector<int> exps;
    exps.reserve(f_pf.size());
    for (auto [p, e] : f_pf) {
        primes.push_back(p);
        exps.push_back(e);
    }

    i128 S = 0;
    const int k = (int)primes.size();
    for (int mask = 0; mask < (1 << k); ++mask) {
        u64 d = 1;
        int mu = 1;
        u128 sigma = 1;
        for (int i = 0; i < k; ++i) {
            int e = exps[i];
            if (mask & (1 << i)) {
                d *= primes[i];
                mu = -mu;
                --e;
            }
            sigma *= sigma_prime_power(primes[i], e);
        }
        int chi = jacobi_symbol(D, (i64)d);
        if (chi == 0) continue;
        S += (i128)mu * (i128)chi * (i128)sigma;
    }

    const i128 t = (i128)12 * (i128)(2 * h) * (i128)one_minus * S;
    assert(t % w == 0);
    const i128 r3 = t / w;
    assert(r3 >= 0);
    return (u64)r3;
}

static i64 brute_r3(i64 N) {
    int lim = (int)std::sqrt((long double)N);
    std::vector<unsigned char> is_sq((size_t)N + 1, 0);
    for (int i = 0; i <= lim; ++i) is_sq[(size_t)i * i] = 1;

    i64 cnt = 0;
    for (int x = -lim; x <= lim; ++x) {
        i64 x2 = (i64)x * x;
        for (int y = -lim; y <= lim; ++y) {
            i64 rem = N - x2 - (i64)y * y;
            if (rem < 0) continue;
            if (!is_sq[(size_t)rem]) continue;
            cnt += (rem == 0) ? 1 : 2;
        }
    }
    return cnt;
}

static u64 G(u64 n) {
    const u64 N = 8 * n + 3;
    const u64 r3 = r3_sum_of_three_squares(N);
    assert(r3 % 8 == 0);
    return r3 / 8;
}

int main() {
    assert(r3_sum_of_three_squares(75) == (u64)brute_r3(75));
    assert(r3_sum_of_three_squares(8003) == (u64)brute_r3(8003));

    assert(G(9) == 7);
    assert(G(1000) == 78);
    assert(G(1000000ULL) == 2106);

    const u64 n = 17526ULL * 1000000000ULL;
    std::cout << G(n) << "\n";
    return 0;
}
