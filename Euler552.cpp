#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;

static std::vector<u32> primes_up_to(u32 n) {
    if (n < 2) return {};
    std::vector<bool> is_comp(static_cast<std::size_t>(n + 1), false);
    std::vector<u32> primes;
    for (u32 i = 2; i <= n; ++i) {
        if (!is_comp[i]) {
            primes.push_back(i);
            if (i <= n / i) {
                for (u32 j = i * i; j <= n; j += i) is_comp[j] = true;
            }
        }
    }
    return primes;
}

static i64 ext_gcd(i64 a, i64 b, i64& x, i64& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i64 x1 = 0, y1 = 0;
    i64 g = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static u32 mod_inv_prime(u32 a, u32 p) {
    // p is prime, a != 0 mod p.
    i64 x = 0, y = 0;
    i64 g = ext_gcd(static_cast<i64>(a), static_cast<i64>(p), x, y);
    assert(g == 1);
    i64 r = x % static_cast<i64>(p);
    if (r < 0) r += static_cast<i64>(p);
    return static_cast<u32>(r);
}

static u64 solve_S(u32 limit) {
    const std::vector<u32> primes = primes_up_to(limit);
    const std::size_t m = primes.size();
    if (m == 0) return 0;

    // Work with all primes <= limit in parallel.
    // At stage i (1-indexed), we have A_i modulo each prime in `primes`,
    // and M_i = p_1...p_i modulo each prime (where p_j are the global primes).
    std::vector<u32> xmod(m, 1U);  // A_1 = 1
    std::vector<u32> mmod(m, 0U);
    for (std::size_t j = 0; j < m; ++j) {
        mmod[j] = static_cast<u32>(2ULL % primes[j]);  // M_1 = 2
    }

    std::vector<std::uint8_t> hit(m, 0);

    // i is 1-indexed index within the CRT system, but we only build it as far as we can with primes<=limit.
    for (std::size_t i = 1; i <= m; ++i) {
        // Check divisibility of current A_i by primes not yet in the system: indices j >= i.
        for (std::size_t j = i; j < m; ++j) {
            if (!hit[j] && xmod[j] == 0) hit[j] = 1;
        }
        if (i == m) break;

        const u32 p = primes[i];  // p_{i+1}
        const u32 rhs = static_cast<u32>(i + 1);  // residue for p_{i+1}
        const u32 xp = xmod[i];
        const u32 mp = mmod[i];
        assert(mp != 0);
        const u32 inv_mp = mod_inv_prime(mp, p);
        const u32 diff = (rhs >= xp) ? (rhs - xp) : (rhs + p - xp);
        const u32 t = static_cast<u32>((static_cast<u64>(diff) * inv_mp) % p);

        // Update A and modulus for all primes >= p_{i+1}.
        for (std::size_t j = i; j < m; ++j) {
            const u32 q = primes[j];
            const u64 add = (static_cast<u64>(mmod[j]) * t) % q;
            u32 nx = xmod[j] + static_cast<u32>(add);
            if (nx >= q) nx -= q;
            xmod[j] = nx;
            mmod[j] = static_cast<u32>((static_cast<u64>(mmod[j]) * p) % q);
        }
        // After adding modulus p, M becomes divisible by p, so mmod[i] must be 0.
        assert(mmod[i] == 0);
        assert(xmod[i] == rhs % p);
    }

    u64 sum = 0;
    for (std::size_t j = 0; j < m; ++j) {
        if (hit[j]) sum += primes[j];
    }
    return sum;
}

}  // namespace

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Validation point from the problem statement.
    assert(solve_S(50) == 69ULL);

    std::cout << solve_S(300'000U) << '\n';
    return 0;
}

