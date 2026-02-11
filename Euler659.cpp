#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod18 = 1'000'000'000'000'000'000ULL;

u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1ULL) result = mod_mul(result, base, mod);
        base = mod_mul(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

u64 tonelli_shanks(u64 n, u64 p) {
    if (n == 0) return 0;
    if (p == 2) return n;
    if (mod_pow(n, (p - 1) / 2, p) != 1) return 0;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);

    u64 q = p - 1;
    int s = 0;
    while ((q & 1ULL) == 0) {
        q >>= 1ULL;
        ++s;
    }

    u64 z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) ++z;

    u64 m = static_cast<u64>(s);
    u64 c = mod_pow(z, q, p);
    u64 t = mod_pow(n, q, p);
    u64 r = mod_pow(n, (q + 1) / 2, p);

    while (t != 1) {
        u64 tt = t;
        u64 i = 0;
        while (tt != 1) {
            tt = mod_mul(tt, tt, p);
            ++i;
        }
        u64 b = mod_pow(c, 1ULL << (m - i - 1), p);
        r = mod_mul(r, b, p);
        c = mod_mul(b, b, p);
        t = mod_mul(t, c, p);
        m = i;
    }

    return r;
}

std::vector<int> sieve_primes(int n) {
    std::vector<int> primes;
    std::vector<bool> composite(static_cast<std::size_t>(n + 1), false);
    for (int i = 2; i <= n; ++i) {
        if (!composite[static_cast<std::size_t>(i)]) primes.push_back(i);
        for (int p : primes) {
            const long long v = 1LL * p * i;
            if (v > n) break;
            composite[static_cast<std::size_t>(v)] = true;
            if (i % p == 0) break;
        }
    }
    return primes;
}

u64 solve_case(int kmax) {
    const int limit = 2 * kmax;
    const std::vector<int> primes = sieve_primes(limit);

    std::vector<u64> rem(static_cast<std::size_t>(kmax + 1), 0);
    std::vector<u32> last_small(static_cast<std::size_t>(kmax + 1), 1);
    for (int k = 1; k <= kmax; ++k) {
        rem[static_cast<std::size_t>(k)] = 4ULL * static_cast<u64>(k) * static_cast<u64>(k) + 1ULL;
    }

    for (int p : primes) {
        if (p == 2 || (p & 3) != 1) continue;

        const u64 root = tonelli_shanks(static_cast<u64>(p - 1), static_cast<u64>(p));
        if (root == 0) continue;

        const u64 inv2 = static_cast<u64>(p + 1) / 2ULL;
        const int k1 = static_cast<int>(mod_mul(root, inv2, static_cast<u64>(p)));
        const int k2 = (k1 == 0) ? 0 : (p - k1);

        auto process_residue = [&](int residue) {
            int start = residue;
            if (start <= 0) start += p;
            for (int k = start; k <= kmax; k += p) {
                u64& v = rem[static_cast<std::size_t>(k)];
                if (v % static_cast<u64>(p) != 0) continue;
                do {
                    v /= static_cast<u64>(p);
                } while (v % static_cast<u64>(p) == 0);
                last_small[static_cast<std::size_t>(k)] = static_cast<u32>(p);
            }
        };

        process_residue(k1);
        if (k2 != k1) process_residue(k2);
    }

    u64 ans = 0;
    for (int k = 1; k <= kmax; ++k) {
        const u64 tail = rem[static_cast<std::size_t>(k)];
        u64 lpf = static_cast<u64>(last_small[static_cast<std::size_t>(k)]);
        if (tail > lpf) lpf = tail;
        ans += lpf % kMod18;
        if (ans >= kMod18) ans %= kMod18;
    }
    return ans % kMod18;
}

}  // namespace

int main() {
    assert(solve_case(1) == 5ULL);
    assert(solve_case(2) == 22ULL);
    assert(solve_case(10) == 1070ULL);
    assert(solve_case(100) == 433752ULL);

    const u64 ans = solve_case(10'000'000);
    std::cout << std::setw(18) << std::setfill('0') << ans << "\n";
    return 0;
}
