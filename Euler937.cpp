#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
constexpr u64 kMod = 1'000'000'007ULL;

std::vector<int> primes_until(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    is_prime[0] = false;
    is_prime[1] = false;
    for (int p = 2; static_cast<u64>(p) * static_cast<u64>(p) <= static_cast<u64>(n); ++p) {
        if (!is_prime[p]) {
            continue;
        }
        for (int q = p * p; q <= n; q += p) {
            is_prime[q] = false;
        }
    }
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(n / 10));
    for (int i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

u64 solve(const int n) {
    const auto primes = primes_until(n);
    std::vector<std::uint8_t> diff(static_cast<std::size_t>(n + 1), 0U);

    for (int p : primes) {
        const int r = p & 7;
        if (r == 1 || r == 3) {
            continue;
        }

        int e = 0;
        int pe = __builtin_parity(static_cast<unsigned>(e));
        for (int q = p; q <= n; q += p) {
            int x = q / p;
            int cnt = 1;
            while (x % p == 0) {
                ++cnt;
                x /= p;
            }
            e += cnt;
            const int ce = __builtin_parity(static_cast<unsigned>(e));
            diff[static_cast<std::size_t>(q)] ^= static_cast<std::uint8_t>(pe ^ ce);
            pe = ce;
        }
    }

    std::uint8_t parity = 0U;
    u64 ans = 0;
    u64 fact = 1;
    for (int k = 1; k <= n; ++k) {
        parity ^= diff[static_cast<std::size_t>(k)];
        fact = (fact * static_cast<u64>(k)) % kMod;
        if (parity == 0U) {
            ans += fact;
            if (ans >= kMod) {
                ans -= kMod;
            }
        }
    }

    return ans;
}

void run_validations() {
    assert(solve(4) == 25ULL);
    assert(solve(7) == 745ULL);
    assert(solve(100) == 709'772'949ULL);
}

}  // namespace

int main() {
    run_validations();
    constexpr int kN = 100'000'000;
    std::cout << solve(kN) << '\n';
    return 0;
}
