#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct OddPrimeTable {
    // index i stores primality of number (2*i + 1)
    std::vector<std::uint8_t> odd_is_prime;

    explicit OddPrimeTable(const int limit) : odd_is_prime(static_cast<std::size_t>((limit + 1) / 2), 1U) {
        if (!odd_is_prime.empty()) {
            odd_is_prime[0] = 0U;  // 1 is not prime
        }
        const int max_i = static_cast<int>(odd_is_prime.size()) - 1;
        for (int i = 1; (2 * i + 1) * (2 * i + 1) <= max_i * 2 + 1; ++i) {
            if (!odd_is_prime[static_cast<std::size_t>(i)]) {
                continue;
            }
            const int p = 2 * i + 1;
            int start = (p * p) >> 1;
            for (int j = start; j <= max_i; j += p) {
                odd_is_prime[static_cast<std::size_t>(j)] = 0U;
            }
        }
    }

    bool is_prime(const int x) const {
        if (x == 2) {
            return true;
        }
        if (x < 2 || (x & 1) == 0) {
            return false;
        }
        return odd_is_prime[static_cast<std::size_t>(x >> 1)] != 0U;
    }
};

u64 solve(const int n) {
    OddPrimeTable prime_table(n);
    const int limit = static_cast<int>(std::sqrt(static_cast<long double>(n)));

    u64 sum = 0ULL;
    for (int v = 2; v <= limit; ++v) {
        const int v2 = v * v;
        const int max_d = n / v2;
        for (int u = 1; u < v; ++u) {
            if (std::gcd(u, v) != 1) {
                continue;
            }

            const int u2 = u * u;
            const int uv = u * v;

            // Odd d can only work when a=2, i.e. d*u^2=3 -> u=1, d=3.
            if (u == 1 && (v % 2 == 0) && max_d >= 3) {
                const int d = 3;
                const int a = u2 * d - 1;
                const int b = uv * d - 1;
                const int c = v2 * d - 1;
                if (prime_table.is_prime(a) && prime_table.is_prime(b) && prime_table.is_prime(c)) {
                    sum += static_cast<u64>(a + b + c);
                }
            }

            for (int d = 2; d <= max_d; d += 2) {
                const int a = u2 * d - 1;
                const int b = uv * d - 1;
                const int c = v2 * d - 1;
                if (prime_table.is_prime(a) && prime_table.is_prime(b) && prime_table.is_prime(c)) {
                    sum += static_cast<u64>(a + b + c);
                }
            }
        }
    }
    return sum;
}

u64 brute(const int n) {
    OddPrimeTable prime_table(n);
    std::vector<int> primes;
    for (int x = 2; x < n; ++x) {
        if (prime_table.is_prime(x)) {
            primes.push_back(x);
        }
    }

    u64 sum = 0ULL;
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const int a = primes[i];
        for (std::size_t j = i + 1; j < primes.size(); ++j) {
            const int b = primes[j];
            const i64 lhs = static_cast<i64>(b + 1) * static_cast<i64>(b + 1);
            for (std::size_t k = j + 1; k < primes.size(); ++k) {
                const int c = primes[k];
                const i64 rhs = static_cast<i64>(a + 1) * static_cast<i64>(c + 1);
                if (lhs == rhs) {
                    sum += static_cast<u64>(a + b + c);
                }
            }
        }
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(100) != 1'035ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    if (solve(500) != brute(500)) {
        std::cerr << "Checkpoint failed: fast/brute mismatch at 500\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr int n = 100'000'000;
    std::cout << solve(n) << '\n';
    return 0;
}
