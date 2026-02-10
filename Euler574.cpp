#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Project Euler 574: Verifying Primes
//
// For a prime p, a verification is given by a prime q and coprime integers A>=B>0 such that:
// - AB is divisible by every prime < q
// - either p = A+B < q^2  (sum) or p = A-B < q^2 (difference)
//
// Let V(p) be the minimum possible A among all verifications, and S(n)=sum_{p<n, p prime} V(p).
//
// Key observations:
// 1) In any verification, AB is divisible by the primorial P = Π_{r<q} r and gcd(A,B)=1,
//    so every prime r<q divides exactly one of A or B.
//    Hence we can write A = D*x, B = E*y where D*E = P and gcd(D,E)=1 (a partition of primes < q).
// 2) For p prime, in the sum case (p=A+B), gcd(A,B)=gcd(A,p)=1 automatically since 0<A<p.
//    In the difference case (p=A-B), gcd(A,B)=gcd(A,p), so gcd(A,B)=1 iff p ∤ A.
// 3) For fixed (q, D, E), both sum and difference impose the same congruences on A:
//      A ≡ 0 (mod D),  A ≡ p (mod E)   (since B is p-A or A-p and must be divisible by E).
//    With D and E coprime, this has a unique solution A0 modulo P.
// 4) If a verification exists for some q, it also exists for any smaller q' with p<q'^2
//    (the divisibility requirement becomes weaker). Therefore the minimal A uses the smallest
//    prime q such that q^2 > p.
//
// So for each prime p we fix q = min{prime : q^2 > p}, build P, enumerate all partitions D|P
// (E=P/D), solve the CRT once per partition, and take the best A among feasible sum/difference
// candidates.

using i128 = __int128_t;
using u64 = std::uint64_t;

static std::string to_string_i128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x > 0) {
        int digit = int(x % 10);
        s.push_back(char('0' + digit));
        x /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

static std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    is_prime.assign(n + 1, true);
    if (n >= 0) is_prime[0] = false;
    if (n >= 1) is_prime[1] = false;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[j] = false;
    }
    std::vector<int> ps;
    for (int i = 2; i <= n; ++i)
        if (is_prime[i]) ps.push_back(i);
    return ps;
}

static int next_prime_after(int x, const std::vector<int>& primes) {
    auto it = std::upper_bound(primes.begin(), primes.end(), x);
    return (it == primes.end()) ? -1 : *it;
}

static i128 egcd(i128 a, i128 b, i128& x, i128& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i128 x1 = 0, y1 = 0;
    i128 g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static i128 mod_inv(i128 a, i128 mod) {
    // mod > 1, gcd(a,mod)=1
    i128 x = 0, y = 0;
    i128 g = egcd(a, mod, x, y);
    (void)g;
    x %= mod;
    if (x < 0) x += mod;
    return x;
}

struct DivInfo {
    i128 D = 1;
    i128 E = 1;
    i128 invD_mod_E = 0;  // only meaningful if E>1
};

static std::vector<DivInfo> build_divisors_with_inverses(const std::vector<int>& primes_lt_q) {
    i128 P = 1;
    for (int r : primes_lt_q) P *= (i128)r;

    std::vector<i128> divisors = {1};
    for (int r : primes_lt_q) {
        const int sz = (int)divisors.size();
        divisors.reserve(sz * 2);
        for (int i = 0; i < sz; ++i) divisors.push_back(divisors[i] * (i128)r);
    }

    std::vector<DivInfo> out;
    out.reserve(divisors.size());
    for (i128 D : divisors) {
        DivInfo di;
        di.D = D;
        di.E = P / D;
        if (di.E > 1) di.invD_mod_E = mod_inv(D % di.E, di.E);
        out.push_back(di);
    }
    return out;
}

static i128 V_of_prime(int p, int q, const std::vector<int>& primes_lt_q, const std::vector<DivInfo>& divs) {
    i128 P = 1;
    for (int r : primes_lt_q) P *= (i128)r;

    i128 best = (i128)1 << 120;  // huge
    const int half = (p + 1) / 2;

    for (const auto& di : divs) {
        const i128 D = di.D;
        const i128 E = di.E;

        // Solve A ≡ 0 (mod D), A ≡ p (mod E) as A0 = D * t, t ≡ p * inv(D) (mod E).
        i128 A0 = 0;
        if (E == 1) {
            A0 = 0;
        } else {
            i128 t = ((i128)p % E) * di.invD_mod_E % E;
            A0 = D * t;  // 0 <= A0 < D*E = P
        }

        // Sum case: p = A + B with A >= B > 0.
        // Solutions are A = A0 + kP. Take the smallest A in [half, p).
        {
            i128 As = A0;
            if (As < half) {
                const i128 need = (i128)half - As;
                const i128 k = (need + P - 1) / P;
                As += k * P;
            }
            if (As < p) best = std::min(best, As);
        }

        // Difference case: p = A - B with B > 0 (i.e. A > p) and gcd(A,B)=1 iff p ∤ A.
        i128 A = A0;
        if (A <= p) {
            i128 k = (p - A) / P + 1;
            A += k * P;
        }
        if (A % p == 0) A += P;
        best = std::min(best, A);
    }

    return best;
}

static i128 S(int n) {
    const auto primes = primes_up_to(n);
    const auto all_primes = primes_up_to(5000);  // enough for sqrt(3800)

    // Cache divisor data per q.
    struct CacheEntry {
        int q = 0;
        std::vector<int> primes_lt_q;
        std::vector<DivInfo> divs;
    };
    std::vector<CacheEntry> cache;

    auto get_cache = [&](int q) -> CacheEntry& {
        for (auto& e : cache)
            if (e.q == q) return e;
        CacheEntry e;
        e.q = q;
        for (int r : all_primes) {
            if (r >= q) break;
            e.primes_lt_q.push_back(r);
        }
        e.divs = build_divisors_with_inverses(e.primes_lt_q);
        cache.push_back(std::move(e));
        return cache.back();
    };

    i128 sum = 0;
    for (int p : primes) {
        if (p >= n) break;
        const int q = next_prime_after((int)std::sqrt((double)p), all_primes);
        auto& e = get_cache(q);
        sum += V_of_prime(p, q, e.primes_lt_q, e.divs);
    }
    return sum;
}

int main() {
    // Validation from the statement.
    if (to_string_i128(S(10)) != "10") {
        std::cerr << "Validation failed: S(10)\n";
        return 1;
    }
    if (to_string_i128(S(200)) != "7177") {
        std::cerr << "Validation failed: S(200)\n";
        return 1;
    }

    std::cout << to_string_i128(S(3800)) << "\n";
    return 0;
}
