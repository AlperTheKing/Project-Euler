#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;

u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / (x + 1ULL)) {
        ++x;
    }
    while (x > 0ULL && x > n / x) {
        --x;
    }
    return x;
}

u64 icbrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::cbrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / ((x + 1ULL) * (x + 1ULL))) {
        ++x;
    }
    while (x > 0ULL && x > n / (x * x)) {
        --x;
    }
    return x;
}

u64 iroot4_u64(u64 n) {
    u64 x = static_cast<u64>(
        std::sqrt(static_cast<long double>(std::sqrt(static_cast<long double>(n)))));
    auto pow4 = [](u64 y) -> u128 { return static_cast<u128>(y) * y * y * y; };
    while (pow4(x + 1ULL) <= n) {
        ++x;
    }
    while (x > 0ULL && pow4(x) > n) {
        --x;
    }
    return x;
}

class PrimeCounting {
public:
    PrimeCounting() { build_sieve(5'000'000); }

    u64 pi(u64 n) {
        if (n <= static_cast<u64>(sieve_limit_)) {
            return pi_small_[static_cast<std::size_t>(n)];
        }
        auto it = pi_cache_.find(n);
        if (it != pi_cache_.end()) {
            return it->second;
        }

        const u64 a = pi(iroot4_u64(n));
        const u64 b = pi(isqrt_u64(n));
        const u64 c = pi(icbrt_u64(n));

        i64 sum = static_cast<i64>(phi(n, static_cast<int>(a))) +
                  static_cast<i64>((b + a - 2ULL) * (b - a + 1ULL) / 2ULL);

        for (u64 i = a + 1ULL; i <= b; ++i) {
            const u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(i - 1ULL)]);
            const u64 w = n / p;
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
        pi_cache_.emplace(n, out);
        return out;
    }

    const std::vector<int>& primes() const { return primes_; }

private:
    void build_sieve(int limit) {
        sieve_limit_ = limit;
        std::vector<bool> is_comp(static_cast<std::size_t>(limit + 1), false);
        pi_small_.assign(static_cast<std::size_t>(limit + 1), 0ULL);

        for (int i = 2; i <= limit; ++i) {
            if (!is_comp[static_cast<std::size_t>(i)]) {
                primes_.push_back(i);
                if (i <= limit / i) {
                    for (int j = i * i; j <= limit; j += i) {
                        is_comp[static_cast<std::size_t>(j)] = true;
                    }
                }
            }
            pi_small_[static_cast<std::size_t>(i)] =
                pi_small_[static_cast<std::size_t>(i - 1)] +
                (!is_comp[static_cast<std::size_t>(i)] ? 1ULL : 0ULL);
        }

        pi_cache_.reserve(1 << 20);
        phi_cache_.reserve(1 << 20);
    }

    u64 phi(u64 x, int s) {
        if (s == 0) {
            return x;
        }
        if (s == 1) {
            return x - x / 2ULL;
        }
        if (s == 2) {
            return x - x / 2ULL - x / 3ULL + x / 6ULL;
        }
        if (s == 3) {
            return x - x / 2ULL - x / 3ULL - x / 5ULL + x / 6ULL + x / 10ULL + x / 15ULL -
                   x / 30ULL;
        }
        if (x <= static_cast<u64>(sieve_limit_) &&
            static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]) >= x) {
            return 1ULL;
        }

        const u64 key = (x << 6U) ^ static_cast<u64>(s);
        auto it = phi_cache_.find(key);
        if (it != phi_cache_.end()) {
            return it->second;
        }

        const u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]);
        const u64 out = phi(x, s - 1) - phi(x / p, s - 1);
        phi_cache_.emplace(key, out);
        return out;
    }

    int sieve_limit_ = 0;
    std::vector<int> primes_;
    std::vector<u64> pi_small_;
    std::unordered_map<u64, u64> pi_cache_;
    std::unordered_map<u64, u64> phi_cache_;
};

u64 count_smooth(u64 n, PrimeCounting& pc) {
    const u64 root = isqrt_u64(n);

    u128 non_smooth = 0;

    for (int p : pc.primes()) {
        if (static_cast<u64>(p) > root) {
            break;
        }
        non_smooth += static_cast<u64>(p);
    }

    const u64 q_max = n / (root + 1ULL);
    for (u64 q = 1; q <= q_max; ++q) {
        const u64 left = n / (q + 1ULL) + 1ULL;
        const u64 right = n / q;
        const u64 l = std::max<u64>(root + 1ULL, left);
        if (l > right) {
            continue;
        }
        const u64 cnt = pc.pi(right) - pc.pi(l - 1ULL);
        non_smooth += static_cast<u128>(q) * static_cast<u128>(cnt);
    }

    return static_cast<u64>(static_cast<u128>(n) - non_smooth);
}

bool is_smooth_bruteforce(u64 n) {
    if (n == 1ULL) {
        return true;
    }
    const long double root = std::sqrt(static_cast<long double>(n));
    u64 x = n;
    for (u64 p = 2; p * p <= x; ++p) {
        if (x % p != 0ULL) {
            continue;
        }
        if (static_cast<long double>(p) >= root) {
            return false;
        }
        while (x % p == 0ULL) {
            x /= p;
        }
    }
    if (x > 1ULL && static_cast<long double>(x) >= root) {
        return false;
    }
    return true;
}

u64 brute_count(u64 n) {
    u64 count = 0;
    for (u64 i = 1; i <= n; ++i) {
        if (is_smooth_bruteforce(i)) {
            ++count;
        }
    }
    return count;
}

}  // namespace

int main() {
    PrimeCounting pc;

    assert(count_smooth(100ULL, pc) == 29ULL);
    assert(count_smooth(1'000ULL, pc) == brute_count(1'000ULL));

    std::cout << count_smooth(10'000'000'000ULL, pc) << "\n";
    return 0;
}
