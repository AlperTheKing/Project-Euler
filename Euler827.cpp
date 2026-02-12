#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using Real = boost::multiprecision::cpp_dec_float_100;

static constexpr u64 kMod = 409120391ULL;

struct Rep {
    bool valid = false;
    Real log2v = 0;
    std::vector<u64> exps;
};

struct BestD {
    bool valid = false;
    Real log2v = 0;
    u64 e2 = 0;
    Rep rep3;
};

static std::mt19937_64 rng(0);

static u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod(r, a, mod);
        }
        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool is_prime64(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return n == p;
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
            if (x == n - 1) {
                comp = false;
                break;
            }
        }
        if (comp) return false;
    }
    return true;
}

static u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2;
    if (n % 3ULL == 0) return 3;
    if (n % 5ULL == 0) return 5;

    std::uniform_int_distribution<u64> dist(2, n - 2);

    while (true) {
        u64 c = dist(rng);
        u64 x = dist(rng);
        u64 y = x;
        u64 d = 1;

        while (d == 1) {
            x = (mul_mod(x, x, n) + c) % n;
            y = (mul_mod(y, y, n) + c) % n;
            y = (mul_mod(y, y, n) + c) % n;
            u64 diff = x > y ? x - y : y - x;
            d = std::gcd(diff, n);
        }

        if (d != n) return d;
    }
}

static void factor_rec(u64 n, std::map<u64, int>& out) {
    if (n == 1) return;
    if (is_prime64(n)) {
        out[n]++;
        return;
    }
    u64 d = pollard_rho(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}

static std::unordered_map<u64, std::vector<u64>> odd_div_cache;

static std::vector<u64> odd_divisors(u64 n) {
    auto it = odd_div_cache.find(n);
    if (it != odd_div_cache.end()) {
        return it->second;
    }

    std::map<u64, int> fac;
    factor_rec(n, fac);

    std::vector<std::pair<u64, int>> vfac(fac.begin(), fac.end());
    std::vector<u64> divs{1};
    for (auto [p, e] : vfac) {
        std::vector<u64> next;
        next.reserve(divs.size() * static_cast<std::size_t>(e + 1));
        u64 pe = 1;
        for (int i = 0; i <= e; ++i) {
            for (u64 d : divs) {
                next.push_back(d * pe);
            }
            pe *= p;
        }
        divs.swap(next);
    }

    std::vector<u64> odd;
    odd.reserve(divs.size());
    for (u64 d : divs) {
        if (d & 1ULL) odd.push_back(d);
    }
    std::sort(odd.begin(), odd.end());
    odd_div_cache[n] = odd;
    return odd;
}

static std::vector<u64> primes_mod4(int residue, int count) {
    std::vector<u64> ps;
    for (u64 x = 2; static_cast<int>(ps.size()) < count; ++x) {
        if (x % 4 != static_cast<u64>(residue)) continue;
        bool prime = true;
        for (u64 p = 2; p * p <= x; ++p) {
            if (x % p == 0) {
                prime = false;
                break;
            }
        }
        if (prime) ps.push_back(x);
    }
    return ps;
}

static const std::vector<u64> P1 = primes_mod4(1, 80);
static const std::vector<u64> P3 = primes_mod4(3, 80);

static std::vector<Real> make_logs(const std::vector<u64>& ps) {
    static const Real ln2 = log(Real(2));
    std::vector<Real> logs(ps.size());
    for (std::size_t i = 0; i < ps.size(); ++i) {
        logs[i] = log(Real(ps[i])) / ln2;
    }
    return logs;
}

static const std::vector<Real> L1 = make_logs(P1);
static const std::vector<Real> L3 = make_logs(P3);
static const Real LOG2_2 = Real(1);

struct DfsKey {
    u64 rem;
    int idx;
    u64 prevf;
    bool operator<(const DfsKey& o) const {
        if (rem != o.rem) return rem < o.rem;
        if (idx != o.idx) return idx < o.idx;
        return prevf < o.prevf;
    }
};

static Rep solve_product(
    u64 P,
    const std::vector<u64>& primes,
    const std::vector<Real>& logs,
    std::unordered_map<u64, Rep>& top_cache
) {
    if (P == 1) {
        Rep r;
        r.valid = true;
        r.log2v = 0;
        return r;
    }

    auto it_top = top_cache.find(P);
    if (it_top != top_cache.end()) {
        return it_top->second;
    }

    std::map<DfsKey, Rep> memo;

    std::function<Rep(u64, int, u64)> dfs = [&](u64 rem, int idx, u64 prevf) -> Rep {
        if (rem == 1) {
            Rep base;
            base.valid = true;
            base.log2v = 0;
            return base;
        }
        if (idx >= static_cast<int>(primes.size())) {
            return Rep{};
        }

        DfsKey key{rem, idx, prevf};
        auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        Rep best;

        auto divs = odd_divisors(rem);
        for (auto rit = divs.rbegin(); rit != divs.rend(); ++rit) {
            u64 f = *rit;
            if (f == 1 || f > prevf) continue;
            u64 e = (f - 1) / 2;
            if (e == 0) continue;

            Rep sub = dfs(rem / f, idx + 1, f);
            if (!sub.valid) continue;

            Rep cand;
            cand.valid = true;
            cand.log2v = logs[idx] * Real(e) + sub.log2v;
            cand.exps.reserve(sub.exps.size() + 1);
            cand.exps.push_back(e);
            cand.exps.insert(cand.exps.end(), sub.exps.begin(), sub.exps.end());

            if (!best.valid || cand.log2v < best.log2v) {
                best = std::move(cand);
            }
        }

        memo[key] = best;
        return best;
    };

    Rep res = dfs(P, 0, P);
    top_cache[P] = res;
    return res;
}

static std::unordered_map<u64, Rep> cache_rep1;
static std::unordered_map<u64, Rep> cache_rep3;
static std::unordered_map<u64, BestD> cache_bestD;

static BestD best_for_D(u64 D) {
    auto it = cache_bestD.find(D);
    if (it != cache_bestD.end()) {
        return it->second;
    }

    BestD best;
    auto divs = odd_divisors(D);
    for (u64 a2 : divs) {
        u64 e2 = 0;
        Real part2 = 0;
        if (a2 > 1) {
            e2 = (a2 + 1) / 2;
            part2 = LOG2_2 * Real(e2);
        }

        u64 C = D / a2;
        Rep rep3 = solve_product(C, P3, L3, cache_rep3);
        if (!rep3.valid) continue;

        Real cur = part2 + rep3.log2v;
        if (!best.valid || cur < best.log2v) {
            best.valid = true;
            best.log2v = cur;
            best.e2 = e2;
            best.rep3 = std::move(rep3);
        }
    }

    cache_bestD[D] = best;
    return best;
}

struct QRep {
    bool valid = false;
    Real log2v = 0;
    Rep rep1;
    u64 e2 = 0;
    Rep rep3;
};

static QRep Q_rep(u64 n) {
    u64 S = n + 1;

    QRep best;
    auto divs = odd_divisors(S);

    for (u64 B : divs) {
        Rep rep1 = solve_product(B, P1, L1, cache_rep1);
        if (!rep1.valid) continue;

        u64 D = static_cast<u64>((static_cast<u128>(2) * S) / B - 1);
        BestD bd = best_for_D(D);
        if (!bd.valid) continue;

        Real cur = rep1.log2v + bd.log2v;
        if (!best.valid || cur < best.log2v) {
            best.valid = true;
            best.log2v = cur;
            best.rep1 = std::move(rep1);
            best.e2 = bd.e2;
            best.rep3 = bd.rep3;
        }
    }

    return best;
}

static u64 rep_mod(const Rep& rep, const std::vector<u64>& primes, u64 mod) {
    u64 r = 1;
    for (std::size_t i = 0; i < rep.exps.size(); ++i) {
        r = mul_mod(r, pow_mod(primes[i], rep.exps[i], mod), mod);
    }
    return r;
}

static u64 Q_mod(u64 n) {
    QRep qr = Q_rep(n);
    u64 r = rep_mod(qr.rep1, P1, kMod);
    r = mul_mod(r, pow_mod(2, qr.e2, kMod), kMod);
    r = mul_mod(r, rep_mod(qr.rep3, P3, kMod), kMod);
    return r;
}

static u128 rep_to_u128(const Rep& rep, const std::vector<u64>& primes) {
    u128 x = 1;
    for (std::size_t i = 0; i < rep.exps.size(); ++i) {
        u64 e = rep.exps[i];
        u128 p = primes[i];
        u128 pw = 1;
        while (e > 0) {
            if (e & 1ULL) pw *= p;
            p *= p;
            e >>= 1ULL;
        }
        x *= pw;
    }
    return x;
}

static u64 Q_exact_small(u64 n) {
    QRep qr = Q_rep(n);
    u128 x = rep_to_u128(qr.rep1, P1);
    x <<= qr.e2;
    x *= rep_to_u128(qr.rep3, P3);
    return static_cast<u64>(x);
}

int main() {
    assert(Q_exact_small(5) == 15ULL);
    assert(Q_exact_small(10) == 48ULL);
    assert(Q_exact_small(1000) == 8064000ULL);

    u64 ans = 0;
    for (int k = 1; k <= 18; ++k) {
        u64 n = 1;
        for (int i = 0; i < k; ++i) n *= 10ULL;
        ans += Q_mod(n);
        ans %= kMod;
    }

    std::cout << ans << '\n';
    return 0;
}
