#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u8 = std::uint8_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'000ULL;

u64 isqrt_u64(const u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1) * (r + 1) <= n) ++r;
    while (r * r > n) --r;
    return r;
}

u64 icbrt_floor_u64(const u64 n) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(n)));
    while ((u128)(r + 1) * (r + 1) * (r + 1) <= n) ++r;
    while ((u128)r * r * r > n) --r;
    return r;
}

struct Sieve {
    int limit = 0;
    std::vector<int> primes;
    std::vector<int> pi;
};

Sieve sieve_with_pi(const int limit) {
    Sieve out;
    out.limit = limit;
    std::vector<u8> is_prime(static_cast<std::size_t>(limit + 1), 1);
    if (limit >= 0) is_prime[0] = 0;
    if (limit >= 1) is_prime[1] = 0;
    for (int p = 2; (static_cast<long long>(p) * p) <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) continue;
        for (int m = p * p; m <= limit; m += p) is_prime[static_cast<std::size_t>(m)] = 0;
    }
    out.pi.assign(static_cast<std::size_t>(limit + 1), 0);
    int cnt = 0;
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            out.primes.push_back(i);
            ++cnt;
        }
        out.pi[static_cast<std::size_t>(i)] = cnt;
    }
    return out;
}

struct Min25PrimeSummatory {
    u64 n = 0;
    u64 sqrt_n = 0;
    std::vector<u64> vals;
    std::vector<u64> g_sum;
    std::vector<int> id_small;
    std::vector<int> id_large;

    int id_of(const u64 x) const {
        if (x <= sqrt_n) return id_small[static_cast<std::size_t>(x)];
        return id_large[static_cast<std::size_t>(n / x)];
    }

    u64 prime_sum_mod(const u64 x) const { return g_sum[static_cast<std::size_t>(id_of(x))]; }
};

Min25PrimeSummatory min25_prime_sum(const u64 n, const std::vector<int>& primes) {
    Min25PrimeSummatory ms;
    ms.n = n;
    ms.sqrt_n = isqrt_u64(n);

    for (u64 l = 1; l <= n;) {
        const u64 v = n / l;
        const u64 r = n / v;
        ms.vals.push_back(v);
        l = r + 1;
    }

    const int m = static_cast<int>(ms.vals.size());
    ms.g_sum.assign(static_cast<std::size_t>(m), 0ULL);
    ms.id_small.assign(static_cast<std::size_t>(ms.sqrt_n + 1), -1);
    ms.id_large.assign(static_cast<std::size_t>(ms.sqrt_n + 1), -1);

    for (int i = 0; i < m; ++i) {
        const u64 v = ms.vals[static_cast<std::size_t>(i)];
        if (v <= ms.sqrt_n) {
            ms.id_small[static_cast<std::size_t>(v)] = i;
        } else {
            ms.id_large[static_cast<std::size_t>(n / v)] = i;
        }
        const u128 vv = static_cast<u128>(v);
        const u128 sum_1_to_v = vv * (vv + 1U) / 2U;
        const u64 sum_mod = static_cast<u64>(sum_1_to_v % static_cast<u128>(kMod));
        ms.g_sum[static_cast<std::size_t>(i)] = (sum_mod + kMod - 1ULL) % kMod;
    }

    u64 primes_sum_mod = 0ULL;
    for (const int p : primes) {
        const u64 pu = static_cast<u64>(p);
        if (pu * pu > n) break;
        for (int i = 0; i < m; ++i) {
            const u64 v = ms.vals[static_cast<std::size_t>(i)];
            if (v < pu * pu) break;
            const u64 vp = v / pu;
            const int j = ms.id_of(vp);
            const u64 sum_j = ms.g_sum[static_cast<std::size_t>(j)];
            const u64 tmp = (sum_j >= primes_sum_mod) ? (sum_j - primes_sum_mod)
                                                      : (sum_j + kMod - primes_sum_mod);
            const u64 delta = (pu % kMod) * tmp % kMod;
            const u64 cur = ms.g_sum[static_cast<std::size_t>(i)];
            ms.g_sum[static_cast<std::size_t>(i)] = (cur >= delta) ? (cur - delta) : (cur + kMod - delta);
        }
        primes_sum_mod += pu % kMod;
        primes_sum_mod %= kMod;
    }
    return ms;
}

struct SmoothCounter {
    const std::vector<int>& primes;
    std::unordered_map<u64, u64> memo;

    explicit SmoothCounter(const std::vector<int>& primes_) : primes(primes_) { memo.reserve(1 << 20); }

    u64 psi(const u64 x, const int idx) {
        if (x == 0) return 0;
        if (idx < 0) return 1;
        if (static_cast<u64>(primes[static_cast<std::size_t>(idx)]) > x) return x;
        const u64 key = (x << 10) | static_cast<u64>(idx);
        const auto it = memo.find(key);
        if (it != memo.end()) return it->second;
        const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(idx)]);
        const u64 res = psi(x, idx - 1) + psi(x / p, idx);
        memo.emplace(key, res);
        return res;
    }
};

u64 psi_easy(const u64 x, const int p, const std::vector<int>& pi) {
    u64 bad = 0;
    u64 l = static_cast<u64>(p) + 1ULL;
    while (l <= x) {
        const u64 q = x / l;
        u64 r = x / q;
        if (r > x) r = x;
        const int cnt = pi[static_cast<std::size_t>(r)] - pi[static_cast<std::size_t>(l - 1)];
        bad += static_cast<u64>(cnt) * q;
        l = r + 1;
    }
    return x - bad;
}

u64 brute(const u32 n) {
    std::vector<u32> lp(static_cast<std::size_t>(n + 1), 0U);
    for (u32 i = 2; i <= n; ++i) {
        if (lp[i] != 0) continue;
        for (u32 j = i; j <= n; j += i) lp[j] = i;
    }
    u64 sum = 0;
    for (u32 i = 2; i <= n; ++i) sum += static_cast<u64>(lp[i]);
    return sum % kMod;
}

u64 solve(const u64 n) {
    if (n < 2) return 0;
    const u64 sqrt_n = isqrt_u64(n);
    const u64 cbrt_floor = icbrt_floor_u64(n);
    u64 cbrt_ceil = cbrt_floor;
    if ((u128)cbrt_ceil * cbrt_ceil * cbrt_ceil < n) ++cbrt_ceil;
    const u64 limit_x = n / cbrt_ceil;
    const u64 sieve_limit_u = std::max(limit_x, sqrt_n);
    const int sieve_limit = static_cast<int>(sieve_limit_u);

    const Sieve sv = sieve_with_pi(sieve_limit);
    const Min25PrimeSummatory ms = min25_prime_sum(n, sv.primes);

    const std::size_t upto_sqrt =
        static_cast<std::size_t>(std::upper_bound(sv.primes.begin(), sv.primes.end(), (int)sqrt_n) - sv.primes.begin());

    int small_k = 0;
    while (small_k < (int)upto_sqrt) {
        const u64 p = static_cast<u64>(sv.primes[static_cast<std::size_t>(small_k)]);
        if ((u128)p * p * p >= n) break;
        ++small_k;
    }

    SmoothCounter smooth(sv.primes);

    u64 ans_mod = 0;
    for (std::size_t i = 0; i < upto_sqrt; ++i) {
        const u64 p = static_cast<u64>(sv.primes[i]);
        const u64 x = n / p;
        const u64 psi = (static_cast<int>(i) < small_k) ? smooth.psi(x, static_cast<int>(i))
                                                        : psi_easy(x, static_cast<int>(p), sv.pi);
        ans_mod += (p % kMod) * (psi % kMod) % kMod;
        ans_mod %= kMod;
    }

    const u64 kmax = n / (sqrt_n + 1);
    for (u64 k = 1; k <= kmax; ++k) {
        const u64 hi = n / k;
        const u64 lo = n / (k + 1);
        const u64 low_bound = (lo < sqrt_n) ? sqrt_n : lo;
        if (hi <= low_bound) continue;
        const u64 sum_hi = ms.prime_sum_mod(hi);
        const u64 sum_lo = ms.prime_sum_mod(low_bound);
        const u64 interval = (sum_hi >= sum_lo) ? (sum_hi - sum_lo) : (sum_hi + kMod - sum_lo);
        ans_mod += (k % kMod) * interval % kMod;
        ans_mod %= kMod;
    }

    return ans_mod;
}

}  // namespace

int main() {
    assert(solve(10) == 32);
    assert(solve(100) == 1915);
    assert(solve(10'000) == 10'118'280);
    assert(solve(1'000'000) == brute(1'000'000));

    std::cout << solve(201'820'182'018ULL) << "\n";
    return 0;
}
