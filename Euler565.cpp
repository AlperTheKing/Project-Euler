#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i64 = std::int64_t;

static std::string to_string_u128(u128 value) {
    if (value == 0) return "0";
    std::string s;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static u128 tri(u64 n) {
    // n*(n+1)/2 in 128-bit.
    return static_cast<u128>(n) * static_cast<u128>(n + 1) / 2;
}

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e) {
        if (e & 1) r = (u128)r * a % mod;
        a = (u128)a * a % mod;
        e >>= 1;
    }
    return r;
}

static i64 egcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0, y1 = 0;
    const i64 g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static u64 mod_inv(u64 a, u64 mod) {
    i64 x = 0, y = 0;
    const i64 g = egcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    assert(g == 1);
    i64 res = x % static_cast<i64>(mod);
    if (res < 0) res += static_cast<i64>(mod);
    return static_cast<u64>(res);
}

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) primes.push_back(i);
    }
    return primes;
}

static u64 order_mod_prime(u64 a, u64 p) {
    // p is prime. Return multiplicative order of a mod p (a != 0 mod p).
    assert(p == 2017);
    const u64 mod = p;
    u64 ord = p - 1;  // 2016 = 2^5 * 3^2 * 7
    const u64 factors[] = {2, 3, 7};
    for (u64 f : factors) {
        while (ord % f == 0 && mod_pow(a, ord / f, mod) == 1) ord /= f;
    }
    return ord;
}

// Sieve primes q <= N with q ≡ -1 (mod p), using a sieve over k where q = p*k - 1.
static std::vector<u64> primes_mod_minus_one(u64 N, u64 p, const std::vector<int>& primes_up_to_sqrtN) {
    assert(p == 2017);
    const u64 k_max = (N + 1) / p;  // p*k - 1 <= N  <=>  k <= (N+1)/p
    std::vector<std::uint8_t> is_comp(k_max + 1, 0);

    for (int r : primes_up_to_sqrtN) {
        if (r == static_cast<int>(p)) continue;
        const u64 rr = static_cast<u64>(r) * static_cast<u64>(r);
        if (rr > N) break;

        const u64 inv = mod_inv(p % static_cast<u64>(r), static_cast<u64>(r));  // k ≡ inv (mod r)
        const u64 k_min = (rr + 1 + p - 1) / p;                                 // ceil((r^2+1)/p)

        u64 k = inv;
        if (k < k_min) {
            const u64 step = static_cast<u64>(r);
            const u64 t = (k_min - k + step - 1) / step;
            k += t * step;
        }
        for (; k <= k_max; k += static_cast<u64>(r)) {
            is_comp[k] = 1;
        }
    }

    std::vector<u64> out;
    out.reserve(2'200'000);
    for (u64 k = 1; k <= k_max; ++k) {
        if (is_comp[k]) continue;
        const u64 q = p * k - 1;
        if (q >= 2 && q <= N) out.push_back(q);
    }
    return out;
}

struct BadPrime {
    u64 q = 0;
    std::vector<int> bad_exps;  // exact exponents e with p|sigma(q^e) within q^e<=N
};

static std::vector<BadPrime> bad_primes_other(u64 N, u64 p, int sqrtN, const std::vector<int>& primes_up_to_sqrtN) {
    std::vector<BadPrime> bad;
    for (int q : primes_up_to_sqrtN) {
        if (q > sqrtN) break;
        if (q == static_cast<int>(p)) continue;
        const u64 a = static_cast<u64>(q) % p;
        if (a == 1 || a == p - 1) continue;  // order 1 is irrelevant here; order 2 handled separately

        const u64 ord = order_mod_prime(a, p);

        // max exponent with q^e <= N
        int emax = 0;
        u64 pw = 1;
        while (pw <= N / static_cast<u64>(q)) {
            pw *= static_cast<u64>(q);
            ++emax;
        }

        std::vector<int> exps;
        for (int m = 1;; ++m) {
            const u64 e = ord * static_cast<u64>(m) - 1;
            if (e > static_cast<u64>(emax)) break;
            exps.push_back(static_cast<int>(e));
        }
        if (exps.empty()) continue;

        bad.push_back(BadPrime{static_cast<u64>(q), std::move(exps)});
    }
    return bad;
}

static u128 sum_exact_bad_exp(u64 N, u64 q, int e) {
    // Sum_{n<=N, v_q(n)=e} n = q^e * sum_{m<=N/q^e, q∤m} m.
    u128 qpow = 1;
    for (int i = 0; i < e; ++i) qpow *= q;
    const u64 M = static_cast<u64>(static_cast<u128>(N) / qpow);
    const u64 Mq = M / q;
    const u128 sum_not_div_q = tri(M) - static_cast<u128>(q) * tri(Mq);
    return qpow * sum_not_div_q;
}

static u128 sum_pair_exact(u64 N, u64 q, int eq, u64 r, int er) {
    // Sum_{n<=N, v_q(n)=eq and v_r(n)=er} n.
    u128 qpow = 1;
    for (int i = 0; i < eq; ++i) qpow *= q;
    u128 rpow = 1;
    for (int i = 0; i < er; ++i) rpow *= r;
    const u128 P = qpow * rpow;
    if (P > N) return 0;
    const u64 M = static_cast<u64>(static_cast<u128>(N) / P);

    const u64 Mq = M / q;
    const u64 Mr = M / r;
    const u64 Mqr = M / (q * r);

    const u128 sum_coprime = tri(M) - static_cast<u128>(q) * tri(Mq) - static_cast<u128>(r) * tri(Mr) +
                             static_cast<u128>(q) * static_cast<u128>(r) * tri(Mqr);
    return P * sum_coprime;
}

static u128 S(u64 N, u64 p) {
    const int sqrtN = static_cast<int>(std::sqrt(static_cast<long double>(N)));
    const auto primes = sieve_primes(sqrtN);

    // Order-2 primes: q ≡ -1 (mod p) with exact exponent 1 (since (min such q)^3 > 1e11).
    const auto p2 = primes_mod_minus_one(N, p, primes);

    // Other bad primes within sqrt(N): exact exponents e = ord(q)*m - 1 within range.
    auto other = bad_primes_other(N, p, sqrtN, primes);

    // Singles
    u128 singles = 0;
    for (u64 q : p2) {
        singles += sum_exact_bad_exp(N, q, 1);
    }
    for (const auto& bp : other) {
        for (int e : bp.bad_exps) singles += sum_exact_bad_exp(N, bp.q, e);
    }

    // Pairs: only need q<r once, but given N=1e11 there are no triple intersections.
    u128 pairs = 0;

    // order2-order2 pairs: q*r <= N, so q must be <= sqrt(N).
    const u64 sqrtN_u = static_cast<u64>(sqrtN);
    const auto it_small_end = std::upper_bound(p2.begin(), p2.end(), sqrtN_u);
    for (auto itq = p2.begin(); itq != it_small_end; ++itq) {
        const u64 q = *itq;
        const u64 lim = N / q;
        auto itr_end = std::upper_bound(itq + 1, p2.end(), lim);
        for (auto itr = itq + 1; itr != itr_end; ++itr) {
            const u64 r = *itr;
            pairs += sum_pair_exact(N, q, 1, r, 1);
        }
    }

    // order2-other pairs: very few due to size constraints (but do it generically).
    for (const auto& bp : other) {
        for (int e : bp.bad_exps) {
            u128 rpow = 1;
            for (int i = 0; i < e; ++i) rpow *= bp.q;
            if (rpow > N) continue;
            const u64 lim = static_cast<u64>(static_cast<u128>(N) / rpow);
            for (u64 q : p2) {
                if (q > lim) break;
                pairs += sum_pair_exact(N, q, 1, bp.q, e);
            }
        }
    }

    // other-other pairs: (empirically none for N=1e11, but keep safe).
    for (std::size_t i = 0; i < other.size(); ++i) {
        for (std::size_t j = i + 1; j < other.size(); ++j) {
            for (int ei : other[i].bad_exps) {
                for (int ej : other[j].bad_exps) {
                    pairs += sum_pair_exact(N, other[i].q, ei, other[j].q, ej);
                }
            }
        }
    }

    return singles - pairs;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const u64 p = 2017;

    // Validation points from the statement.
    assert(to_string_u128(S(1'000'000ULL, p)) == "150850429");
    assert(to_string_u128(S(1'000'000'000ULL, p)) == "249652238344557");

    std::cout << to_string_u128(S(100'000'000'000ULL, p)) << '\n';
    return 0;
}
