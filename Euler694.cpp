#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

std::vector<int> prime_sieve(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n) + 1U, true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }

    for (int p = 2; static_cast<u64>(p) * static_cast<u64>(p) <= static_cast<u64>(n); ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (u64 q = static_cast<u64>(p) * static_cast<u64>(p); q <= static_cast<u64>(n);
             q += static_cast<u64>(p)) {
            is_prime[static_cast<std::size_t>(q)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

u64 icbrt(const u64 n) {
    u64 lo = 0;
    u64 hi = 2'000'000;
    while (lo + 1 < hi) {
        const u64 mid = lo + (hi - lo) / 2;
        const u128 cube = static_cast<u128>(mid) * static_cast<u128>(mid) * static_cast<u128>(mid);
        if (cube <= static_cast<u128>(n)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

u64 summatory_cube_full_divisors(const u64 n) {
    const u64 lim = icbrt(n);
    const std::vector<int> primes = prime_sieve(static_cast<int>(lim));

    u64 total = 0;

    const auto dfs = [&](auto&& self, const std::size_t idx, const u64 cur) -> void {
        total += n / cur;

        for (std::size_t j = idx; j < primes.size(); ++j) {
            const u64 p = static_cast<u64>(primes[j]);
            const u128 p3 = static_cast<u128>(p) * static_cast<u128>(p) * static_cast<u128>(p);
            if (static_cast<u128>(cur) * p3 > static_cast<u128>(n)) {
                break;
            }

            u64 value = static_cast<u64>(static_cast<u128>(cur) * p3);
            while (value <= n) {
                self(self, j + 1, value);
                if (value > n / p) {
                    break;
                }
                value *= p;
            }
        }
    };

    dfs(dfs, 0U, 1U);
    return total;
}

bool is_cube_full(u64 x) {
    if (x == 1U) {
        return true;
    }

    for (u64 p = 2; p * p <= x; ++p) {
        if (x % p != 0U) {
            continue;
        }
        int e = 0;
        while (x % p == 0U) {
            x /= p;
            ++e;
        }
        if (e < 3) {
            return false;
        }
    }

    if (x > 1U) {
        return false;
    }
    return true;
}

u64 brute_s(const u64 n) {
    u64 out = 0;
    for (u64 d = 1; d <= n; ++d) {
        if (n % d == 0U && is_cube_full(d)) {
            ++out;
        }
    }
    return out;
}

u64 brute_S(const u64 n) {
    u64 out = 0;
    for (u64 i = 1; i <= n; ++i) {
        out += brute_s(i);
    }
    return out;
}

}  // namespace

int main() {
    assert(summatory_cube_full_divisors(16U) == 19U);
    assert(summatory_cube_full_divisors(100U) == 126U);
    assert(summatory_cube_full_divisors(10'000U) == 13'344U);
    assert(summatory_cube_full_divisors(250U) == brute_S(250U));

    std::cout << summatory_cube_full_divisors(1'000'000'000'000'000'000ULL) << '\n';
    return 0;
}

