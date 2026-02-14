#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <functional>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr u64 kDenTarget = 20010ULL;
static constexpr u64 kL = 308ULL;  // lcm(4,22,28)

static std::vector<int> sieve_primes(int n) {
    std::vector<bool> is_comp(static_cast<std::size_t>(n + 1), false);
    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (!is_comp[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
            if (i <= n / i) {
                for (int j = i * i; j <= n; j += i) is_comp[static_cast<std::size_t>(j)] = true;
            }
        }
    }
    return primes;
}

static bool is_prime_u64(u64 x, const std::vector<int>& primes) {
    if (x < 2) return false;
    for (int p : primes) {
        const u64 pp = static_cast<u64>(p);
        if (pp * pp > x) break;
        if (x % pp == 0) return x == pp;
    }
    return true;
}

static int mod_inv(int a, int mod) {
    // Extended Euclid, mod is prime in our use.
    int t = 0, newt = 1;
    int r = mod, newr = a;
    while (newr != 0) {
        const int q = r / newr;
        const int tmp_t = t - q * newt;
        t = newt;
        newt = tmp_t;
        const int tmp_r = r - q * newr;
        r = newr;
        newr = tmp_r;
    }
    if (r != 1) return 0;
    if (t < 0) t += mod;
    return t;
}

static u64 bernoulli_denominator_even_k(u64 k, const std::vector<int>& primes_small) {
    // Von Staudt–Clausen: denom(B_k) = Π_{p-1 | k} p for even k>=2.
    // Compute by enumerating divisors d|k and checking if d+1 is prime.
    std::vector<std::pair<u64, int>> fac;
    u64 x = k;
    for (int p : primes_small) {
        const u64 pp = static_cast<u64>(p);
        if (pp * pp > x) break;
        if (x % pp == 0) {
            int e = 0;
            while (x % pp == 0) {
                x /= pp;
                ++e;
            }
            fac.emplace_back(pp, e);
        }
    }
    if (x > 1) fac.emplace_back(x, 1);

    u64 den = 1;
    std::vector<u64> divs{1};
    for (const auto& [p, e] : fac) {
        const std::size_t cur = divs.size();
        u64 pe = 1;
        for (int i = 1; i <= e; ++i) {
            pe *= p;
            for (std::size_t j = 0; j < cur; ++j) divs.push_back(divs[j] * pe);
        }
    }
    for (u64 d : divs) {
        const u64 cand = d + 1;
        if (is_prime_u64(cand, primes_small)) den *= cand;
    }
    return den;
}

static u64 brute_F_10() {
    // Brute the first 10 k with D(k)=20010 by scanning multiples of 308.
    const int lim = 200000;  // enough to reach F(10)=96404.
    const auto primes = sieve_primes(static_cast<int>(std::sqrt(lim + 1)) + 5);

    int found = 0;
    u64 tenth = 0;
    for (u64 m = 1; kL * m <= static_cast<u64>(lim); ++m) {
        const u64 k = kL * m;
        const u64 den = bernoulli_denominator_even_k(k, primes);
        if (den == kDenTarget) {
            ++found;
            if (found == 10) {
                tenth = k;
                break;
            }
        }
    }
    return tenth;
}

static u64 kth_multiplier_with_den(u64 M, u64 kth) {
    // m is valid iff no prime p outside {2,3,5,23,29} has (p-1) | (308*m).
    // For each such p, let r=(p-1)/gcd(p-1,308). Then r | m implies invalid.
    //
    // We compute all forbidden r<=M by sieving primes of the form p=g*u+1 with g|308.

    const u64 pmax = kL * M + 1ULL;
    const int qmax = static_cast<int>(std::sqrt(static_cast<long double>(pmax))) + 2;
    const std::vector<int> primes = sieve_primes(qmax);

    std::vector<char> forbidden(static_cast<std::size_t>(M + 1), 0);

    auto is_target_prime = [](u64 p) -> bool {
        return p == 2ULL || p == 3ULL || p == 5ULL || p == 23ULL || p == 29ULL;
    };

    const int gs[] = {2, 4, 14, 22, 28, 44, 154, 308};
    for (int g : gs) {
        const u64 maxp_g = static_cast<u64>(g) * M + 1ULL;
        const u64 qlim = static_cast<u64>(std::sqrt(static_cast<long double>(maxp_g))) + 1ULL;

        std::vector<char> is_prime_u(static_cast<std::size_t>(M + 1), 1);
        is_prime_u[0] = 0;

        for (int q : primes) {
            if (static_cast<u64>(q) > qlim) break;
            if (g % q == 0) continue;

            const int a = g % q;
            const int inv = mod_inv(a, q);
            // u ≡ (-1) * g^{-1} (mod q)
            const int u0 = static_cast<int>((static_cast<long long>(q - 1) * inv) % q);

            const u64 qq = static_cast<u64>(q);
            const u64 u_start = (qq * qq - 1ULL + static_cast<u64>(g) - 1ULL) / static_cast<u64>(g);
            u64 u = static_cast<u64>(u0);
            if (u < u_start) {
                u += ((u_start - u + qq - 1ULL) / qq) * qq;
            }
            for (; u <= M; u += qq) is_prime_u[static_cast<std::size_t>(u)] = 0;
        }

        const u64 div = kL / static_cast<u64>(g);
        for (u64 u = 1; u <= M; ++u) {
            if (!is_prime_u[static_cast<std::size_t>(u)]) continue;
            const u64 p = static_cast<u64>(g) * u + 1ULL;
            if (is_target_prime(p)) continue;
            const u64 r = u / std::gcd(u, div);
            forbidden[static_cast<std::size_t>(r)] = 1;
        }
    }

    std::vector<char> invalid(static_cast<std::size_t>(M + 1), 0);
    for (u64 r = 2; r <= M; ++r) {
        if (!forbidden[static_cast<std::size_t>(r)]) continue;
        for (u64 x = r; x <= M; x += r) invalid[static_cast<std::size_t>(x)] = 1;
    }

    u64 cnt = 0;
    for (u64 m = 1; m <= M; ++m) {
        if (invalid[static_cast<std::size_t>(m)]) continue;
        ++cnt;
        if (cnt == kth) return m;
    }
    return 0ULL;
}

static u64 solve_F_1e5() {
    const u64 kth = 100'000ULL;
    u64 M = 4'000'000ULL;
    while (true) {
        const u64 m = kth_multiplier_with_den(M, kth);
        if (m != 0) return kL * m;
        M *= 2;
    }
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Statement validations.
    {
        const auto primes = sieve_primes(1000);
        assert(bernoulli_denominator_even_k(4, primes) == 30ULL);
        assert(bernoulli_denominator_even_k(308, primes) == kDenTarget);
    }
    assert(brute_F_10() == 96404ULL);

    std::cout << solve_F_1e5() << '\n';
    return 0;
}

