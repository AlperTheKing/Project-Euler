#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1'000'000'007LL;

class PrimeCounter {
  public:
    PrimeCounter() { sieve(); }

    u64 count_primes(const u64 x) { return lehmer(x); }

  private:
    static constexpr int kSieveLimit = 5'000'000;

    std::vector<int> primes_;
    std::vector<int> pi_;
    std::unordered_map<u64, u64> lehmer_cache_;

    void sieve() {
        std::vector<bool> is_prime(static_cast<std::size_t>(kSieveLimit + 1), true);
        is_prime[0] = false;
        is_prime[1] = false;

        for (int i = 2; static_cast<i64>(i) * i <= kSieveLimit; ++i) {
            if (!is_prime[static_cast<std::size_t>(i)]) {
                continue;
            }
            for (int j = i * i; j <= kSieveLimit; j += i) {
                is_prime[static_cast<std::size_t>(j)] = false;
            }
        }

        pi_.assign(static_cast<std::size_t>(kSieveLimit + 1), 0);
        for (int i = 2; i <= kSieveLimit; ++i) {
            if (is_prime[static_cast<std::size_t>(i)]) {
                primes_.push_back(i);
            }
            pi_[static_cast<std::size_t>(i)] =
                pi_[static_cast<std::size_t>(i - 1)] + (is_prime[static_cast<std::size_t>(i)] ? 1 : 0);
        }
    }

    static u64 isqrt(u64 x) {
        u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
        while ((r + 1) <= x / (r + 1)) {
            ++r;
        }
        while (r > x / r) {
            --r;
        }
        return r;
    }

    static u64 icbrt(u64 x) {
        u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
        while ((r + 1) <= x / (r + 1) / (r + 1)) {
            ++r;
        }
        while (r > 0 && r > x / r / r) {
            --r;
        }
        return r;
    }

    static u64 iroot4(u64 x) {
        u64 r = static_cast<u64>(std::sqrt(std::sqrt(static_cast<long double>(x))));
        while ((r + 1) <= x / (r + 1) / (r + 1) / (r + 1)) {
            ++r;
        }
        while (r > 0 && r > x / r / r / r) {
            --r;
        }
        return r;
    }

    u64 phi(const u64 x, const int s) {
        if (s == 0) {
            return x;
        }
        if (s == 1) {
            return x - x / 2;
        }
        if (s == 2) {
            return x - x / 2 - x / 3 + x / 6;
        }
        if (s == 3) {
            return x - x / 2 - x / 3 - x / 5 + x / 6 + x / 10 + x / 15 - x / 30;
        }
        if (x <= kSieveLimit && static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]) *
                                     static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]) >
                                     x) {
            return static_cast<u64>(pi_[static_cast<std::size_t>(x)] - s + 1);
        }
        return phi(x, s - 1) - phi(x / static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]), s - 1);
    }

    u64 lehmer(const u64 x) {
        if (x <= kSieveLimit) {
            return static_cast<u64>(pi_[static_cast<std::size_t>(x)]);
        }

        const auto cached = lehmer_cache_.find(x);
        if (cached != lehmer_cache_.end()) {
            return cached->second;
        }

        const u64 a = lehmer(iroot4(x));
        const u64 b = lehmer(isqrt(x));
        const u64 c = lehmer(icbrt(x));

        u64 sum = phi(x, static_cast<int>(a));
        sum += (b + a - 2) * (b - a + 1) / 2;

        for (u64 i = a + 1; i <= b; ++i) {
            const u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(i - 1)]);
            const u64 w = x / p;
            sum -= lehmer(w);
            if (i <= c) {
                const u64 bi = lehmer(isqrt(w));
                for (u64 j = i; j <= bi; ++j) {
                    const u64 pj = static_cast<u64>(primes_[static_cast<std::size_t>(j - 1)]);
                    sum -= lehmer(w / pj) - (j - 1);
                }
            }
        }

        lehmer_cache_.emplace(x, sum);
        return sum;
    }
};

i64 solve(const u64 n) {
    PrimeCounter prime_counter;

    const u64 y = std::min<u64>(n, 100'000'000ULL);
    const i64 n_mod = static_cast<i64>(n % static_cast<u64>(kMod));
    i64 ans = 0;

    const auto add_q = [&](const u64 q, const u64 count, i64& total) {
        const i64 q_mod = static_cast<i64>(q % static_cast<u64>(kMod));
        i64 term = static_cast<i64>(((__int128)n_mod * q_mod - (__int128)q_mod * q_mod) % kMod);
        if (term < 0) {
            term += kMod;
        }
        const i64 count_mod = static_cast<i64>(count % static_cast<u64>(kMod));
        total = static_cast<i64>((total + (__int128)2 * term % kMod * count_mod) % kMod);
    };

    {
        const u64 segment_size = 1'000'000ULL;
        const int root = static_cast<int>(std::sqrt(static_cast<long double>(y)));

        std::vector<bool> base_mark(static_cast<std::size_t>(root + 1), true);
        std::vector<int> base_primes;
        for (int i = 2; i <= root; ++i) {
            if (!base_mark[static_cast<std::size_t>(i)]) {
                continue;
            }
            base_primes.push_back(i);
            if (static_cast<i64>(i) * i <= root) {
                for (int j = i * i; j <= root; j += i) {
                    base_mark[static_cast<std::size_t>(j)] = false;
                }
            }
        }

        for (u64 low = 2; low <= y; low += segment_size) {
            const u64 high = std::min(y, low + segment_size - 1);
            std::vector<bool> is_prime(static_cast<std::size_t>(high - low + 1), true);

            for (const int p : base_primes) {
                const u64 prime = static_cast<u64>(p);
                u64 start = (low + prime - 1) / prime * prime;
                const u64 p2 = prime * prime;
                if (start < p2) {
                    start = p2;
                }
                if (start > high) {
                    continue;
                }
                for (u64 v = start; v <= high; v += prime) {
                    is_prime[static_cast<std::size_t>(v - low)] = false;
                }
            }

            for (u64 p = low; p <= high; ++p) {
                if (!is_prime[static_cast<std::size_t>(p - low)] || p < 2) {
                    continue;
                }
                u64 power = p;
                while (power <= n) {
                    add_q(n / power, 1, ans);
                    if (power > n / p) {
                        break;
                    }
                    power *= p;
                }
            }
        }
    }

    if (y < n) {
        const u64 q_max = n / (y + 1);
        for (u64 q = 1; q <= q_max; ++q) {
            u64 l = n / (q + 1) + 1;
            u64 r = n / q;
            if (r <= y) {
                continue;
            }
            if (l <= y) {
                l = y + 1;
            }
            if (l > r) {
                continue;
            }

            const u64 count = prime_counter.count_primes(r) - prime_counter.count_primes(l - 1);
            if (count != 0) {
                add_q(q, count, ans);
            }
        }
    }

    return ans;
}

}  // namespace

int main() {
    assert(solve(10) == 210);
    assert(solve(100) == 37'018);

    std::cout << solve(1'000'000'000'000ULL) << '\n';
    return 0;
}
