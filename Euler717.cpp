#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i128 = __int128_t;
using u64 = std::uint64_t;

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>(static_cast<i128>(a) * static_cast<i128>(b) % static_cast<i128>(mod));
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    if (mod == 1) {
        return 0;
    }
    u64 result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

u64 g_value(const u64 p) {
    const u64 e = pow_mod(2, p, p - 1);
    const u64 r = pow_mod(2, e, p);

    const u64 t = (r * ((p + 1) / 2)) % p;

    const u64 mod = p * p;
    const u64 two_p_mod = pow_mod(2, p, mod);
    u64 u = mul_mod(t, two_p_mod, mod);
    if (u < r) {
        u += mod;
    }
    u -= r;

    assert(u % p == 0);
    return u / p;
}

u64 G(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n), true);
    if (n > 0) {
        is_prime[0] = false;
    }
    if (n > 1) {
        is_prime[1] = false;
    }

    for (int i = 2; static_cast<long long>(i) * i < n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j < n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    u64 sum = 0;
    for (int p = 3; p < n; p += 2) {
        if (is_prime[static_cast<std::size_t>(p)]) {
            sum += g_value(static_cast<u64>(p));
        }
    }
    return sum;
}

}  // namespace

int main() {
    assert(g_value(31) == 17);
    assert(G(100) == 474);
    assert(G(10'000) == 2'819'236);

    std::cout << G(10'000'000) << '\n';
    return 0;
}
