#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / (x + 1ULL)) ++x;
    while (x > 0ULL && x > n / x) --x;
    return x;
}

static u64 icbrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::cbrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / ((x + 1ULL) * (x + 1ULL))) ++x;
    while (x > 0ULL && x > n / (x * x)) --x;
    return x;
}

static u64 iroot4_u64(u64 n) {
    const long double x = static_cast<long double>(n);
    u64 r = static_cast<u64>(std::sqrt(std::sqrt(x)));
    auto pow4 = [](u64 y) -> u128 { return static_cast<u128>(y) * y * y * y; };
    while (pow4(r + 1ULL) <= n) ++r;
    while (r > 0ULL && pow4(r) > n) --r;
    return r;
}

class PrimeCounter {
public:
    explicit PrimeCounter(int sieve_limit = 5'000'000) { build_sieve(sieve_limit); }

    u64 pi(u64 n) {
        if (n <= static_cast<u64>(sieve_limit_)) {
            return pi_small_[static_cast<std::size_t>(n)];
        }
        auto it = pi_cache_.find(n);
        if (it != pi_cache_.end()) return it->second;

        const u64 a = pi(iroot4_u64(n));
        const u64 b = pi(isqrt_u64(n));
        const u64 c = pi(icbrt_u64(n));

        i64 sum = static_cast<i64>(phi(n, static_cast<int>(a))) +
                  static_cast<i64>((b + a - 2ULL) * (b - a + 1ULL) / 2ULL);

        for (u64 i = a + 1ULL; i <= b; ++i) {
            const u64 w = n / static_cast<u64>(primes_[static_cast<std::size_t>(i - 1ULL)]);
            sum -= static_cast<i64>(pi(w));
            if (i <= c) {
                const u64 lim = pi(isqrt_u64(w));
                for (u64 j = i; j <= lim; ++j) {
                    const u64 pj = static_cast<u64>(primes_[static_cast<std::size_t>(j - 1ULL)]);
                    sum -= static_cast<i64>(pi(w / pj) - (j - 1ULL));
                }
            }
        }

        const u64 out = static_cast<u64>(sum);
        pi_cache_[n] = out;
        return out;
    }

private:
    void build_sieve(int limit) {
        sieve_limit_ = limit;
        std::vector<bool> is_comp(static_cast<std::size_t>(limit + 1), false);
        pi_small_.assign(static_cast<std::size_t>(limit + 1), 0ULL);
        for (int i = 2; i <= limit; ++i) {
            if (!is_comp[static_cast<std::size_t>(i)]) {
                primes_.push_back(i);
                if (i <= limit / i) {
                    for (int j = i * i; j <= limit; j += i) is_comp[static_cast<std::size_t>(j)] = true;
                }
            }
            pi_small_[static_cast<std::size_t>(i)] =
                pi_small_[static_cast<std::size_t>(i - 1)] + (!is_comp[static_cast<std::size_t>(i)] ? 1ULL : 0ULL);
        }
    }

    u64 phi(u64 x, int s) {
        if (s == 0) return x;
        if (s == 1) return x - x / 2ULL;
        if (s == 2) return x - x / 2ULL - x / 3ULL + x / 6ULL;
        if (s == 3) return x - x / 2ULL - x / 3ULL - x / 5ULL + x / 6ULL + x / 10ULL + x / 15ULL - x / 30ULL;
        if (x <= static_cast<u64>(sieve_limit_) && static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]) >= x) return 1ULL;

        // (x,s) packing: safe here because s <= pi(iroot4(n)) <= 64 for our n-range.
        const u64 key = (x << 6U) ^ static_cast<u64>(s);
        auto it = phi_cache_.find(key);
        if (it != phi_cache_.end()) return it->second;

        const u64 ps = static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]);
        const u64 out = phi(x, s - 1) - phi(x / ps, s - 1);
        phi_cache_[key] = out;
        return out;
    }

    int sieve_limit_ = 0;
    std::vector<int> primes_;
    std::vector<u64> pi_small_;
    std::unordered_map<u64, u64> pi_cache_;
    std::unordered_map<u64, u64> phi_cache_;
};

static u128 S_of(u64 n, PrimeCounter& pc) {
    const u64 pi_n = pc.pi(n);
    const u64 pi_n_minus2 = (n >= 2 ? pc.pi(n - 2) : 0ULL);

    const u64 m = n / 2ULL;
    const u64 even_2prime = (m >= 2 ? (m - 1ULL) : 0ULL);  // even i>=4
    const u64 odd_2prime = (pi_n_minus2 >= 1 ? (pi_n_minus2 - 1ULL) : 0ULL);  // odd i=2+p, p odd prime

    const u128 sum_even_kge3 = (m >= 3 ? (static_cast<u128>(m - 1ULL) * (m - 2ULL)) / 2U : 0U);
    const u64 m2 = (n >= 1 ? (n - 1ULL) / 2ULL : 0ULL);
    const u128 sum_odd_kge3 = (m2 >= 3 ? (static_cast<u128>(m2 - 1ULL) * (m2 - 2ULL)) / 2U : 0U);

    return static_cast<u128>(pi_n) + static_cast<u128>(even_2prime) + static_cast<u128>(odd_2prime) + sum_even_kge3 +
           sum_odd_kge3;
}

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

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    PrimeCounter pc;

    // Validation points from the problem statement.
    assert(S_of(10, pc) == 20U);
    assert(S_of(100, pc) == 2402U);
    assert(S_of(1000, pc) == 248838U);

    u64 f0 = 0, f1 = 1;
    u128 total = 0;
    for (int k = 2; k <= 44; ++k) {
        const u64 f = f0 + f1;
        f0 = f1;
        f1 = f;
        if (k >= 3) total += S_of(f, pc);
    }

    std::cout << to_string_u128(total) << '\n';
    return 0;
}

