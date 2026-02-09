#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

i64 extended_gcd(const i64 a, const i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0;
    i64 y1 = 0;
    const i64 g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

u64 mod_inverse(u64 a, u64 mod) {
    a %= mod;
    if (a == 0ULL) {
        return 0ULL;
    }
    i64 x = 0;
    i64 y = 0;
    const i64 g = extended_gcd(static_cast<i64>(a), static_cast<i64>(mod), x, y);
    if (g != 1) {
        return 0ULL;
    }
    i64 res = x % static_cast<i64>(mod);
    if (res < 0) {
        res += static_cast<i64>(mod);
    }
    return static_cast<u64>(res);
}

std::vector<int> small_primes_up_to(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; p * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int m = p * p; m <= n; m += p) {
            is_prime[static_cast<std::size_t>(m)] = false;
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

std::vector<u64> segmented_primes(const u64 low, const u64 high_exclusive) {
    if (high_exclusive <= low) {
        return {};
    }

    const u64 high_inclusive = high_exclusive - 1ULL;
    int root = 1;
    while (static_cast<u64>(root) * static_cast<u64>(root) <= high_inclusive) {
        ++root;
    }
    --root;

    const std::vector<int> base_primes = small_primes_up_to(root);
    const u64 size = high_exclusive - low;
    std::vector<bool> is_prime(static_cast<std::size_t>(size), true);

    for (const int p_int : base_primes) {
        const u64 p = static_cast<u64>(p_int);
        u64 start = (low + p - 1ULL) / p * p;
        const u64 pp = p * p;
        if (start < pp) {
            start = pp;
        }
        for (u64 x = start; x < high_exclusive; x += p) {
            is_prime[static_cast<std::size_t>(x - low)] = false;
        }
    }

    if (low == 0ULL) {
        is_prime[0] = false;
        if (size > 1ULL) {
            is_prime[1] = false;
        }
    } else if (low == 1ULL) {
        is_prime[0] = false;
    }

    std::vector<u64> primes;
    for (u64 i = 0ULL; i < size; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(low + i);
        }
    }
    return primes;
}

u64 D(const u64 a, const u64 b, const u64 k) {
    const u64 m = k - 1ULL;
    const std::vector<u64> primes = segmented_primes(a, a + b);

    u64 total = 0ULL;
    for (const u64 p : primes) {
        const u64 term = mod_inverse(m, p);  // equals d(p,p-1,k) mod p
        total += term;
    }
    return total;
}

bool run_checkpoints() {
    if (D(101ULL, 1ULL, 10ULL) != 45ULL) {
        std::cerr << "Checkpoint failed: D(101,1,10)\n";
        return false;
    }
    if (D(1'000ULL, 100ULL, 100ULL) != 8'334ULL) {
        std::cerr << "Checkpoint failed: D(10^3,10^2,10^2)\n";
        return false;
    }
    if (D(1'000'000ULL, 1'000ULL, 1'000ULL) != 38'162'302ULL) {
        std::cerr << "Checkpoint failed: D(10^6,10^3,10^3)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 a = 1'000'000'000ULL;
    constexpr u64 b = 100'000ULL;
    constexpr u64 k = 100'000ULL;
    std::cout << D(a, b, k) << '\n';
    return 0;
}
