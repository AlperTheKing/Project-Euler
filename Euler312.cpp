#include <cstdint>
#include <iostream>
#include <string>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }
    base %= mod;
    u64 result = 1ULL % mod;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base, mod);
        }
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }
    return result;
}

u64 pow_u64(u64 base, int exp) {
    u64 result = 1ULL;
    while (exp-- > 0) {
        result *= base;
    }
    return result;
}

u64 inverse_mod(u64 a, const u64 mod) {
    a %= mod;
    // mod is prime-power of 13 in this problem, so a^(phi-1) gives inverse for gcd(a,mod)=1.
    const u64 phi = mod - mod / 13ULL;
    return pow_mod(a, phi - 1ULL, mod);
}

u64 crt_coprime(const u64 a1, const u64 m1, const u64 a2, const u64 m2) {
    const u64 inv_m1 = inverse_mod(m1 % m2, m2);
    const u64 delta = (a2 >= a1 % m2) ? (a2 - (a1 % m2)) : (a2 + m2 - (a1 % m2));
    const u64 t = mul_mod(delta, inv_m1, m2);
    return a1 + m1 * t;
}

u64 c_mod_iterative(const u64 n, const u64 mod) {
    if (mod == 1ULL) {
        return 0ULL;
    }
    if (n <= 2ULL) {
        return 1ULL % mod;
    }

    u64 c = 8ULL % mod;  // C(3)
    for (u64 i = 4ULL; i <= n; ++i) {
        const u64 v = (3ULL * c) % mod;
        c = mul_mod(mul_mod(v, v, mod), v, mod);
    }
    return c;
}

u64 c_mod_13_power_from_n_residue(const u64 n_residue, const bool n_is_large, const int k) {
    const u64 mod = pow_u64(13ULL, k);

    if (!n_is_large && n_residue <= 2ULL) {
        return 1ULL % mod;
    }

    const u64 ord12 = 2ULL * pow_u64(13ULL, k - 1);        // ord_{13^k}(12)
    const u64 two_ord = 2ULL * ord12;                      // 4 * 13^{k-1}
    const u64 lambda = (k >= 2) ? (12ULL * pow_u64(13ULL, k - 2)) : 2ULL;

    const u64 exp = (n_residue + lambda - 2ULL) % lambda;  // n - 2 mod lambda(two_ord)
    const u64 x = pow_mod(3ULL, exp, two_ord);              // x = 3^{n-2} mod (2*ord)

    const u64 numer = (x + two_ord - 3ULL) % two_ord;       // x - 3 is even in this modulus
    const u64 e = (numer / 2ULL) % ord12;                   // ((3^{n-2}-3)/2) mod ord

    return mul_mod(8ULL % mod, pow_mod(12ULL, e, mod), mod);
}

u64 solve() {
    const u64 n0 = 10000ULL;

    // a = C(10000). For n >= 4, C(n) is divisible by 12.
    const u64 a_mod_13_4 = c_mod_13_power_from_n_residue(n0, false, 4);
    const u64 a_mod_12x13_4 = crt_coprime(0ULL, 12ULL, a_mod_13_4, pow_u64(13ULL, 4));

    // b = C(a). Need b mod (12 * 13^6) to evaluate C(b) mod 13^8.
    const u64 b_mod_13_6 = c_mod_13_power_from_n_residue(a_mod_12x13_4, true, 6);
    const u64 b_mod_12x13_6 = crt_coprime(0ULL, 12ULL, b_mod_13_6, pow_u64(13ULL, 6));

    // c = C(b) mod 13^8.
    return c_mod_13_power_from_n_residue(b_mod_12x13_6, true, 8);
}

bool run_checkpoints() {
    // Given checks from statement.
    if (c_mod_iterative(1ULL, 1000000000000000ULL) != 1ULL) {
        std::cerr << "Checkpoint failed: C(1)\n";
        return false;
    }
    if (c_mod_iterative(2ULL, 1000000000000000ULL) != 1ULL) {
        std::cerr << "Checkpoint failed: C(2)\n";
        return false;
    }
    if (c_mod_iterative(5ULL, 1000000000000000ULL) != 71328803586048ULL) {
        std::cerr << "Checkpoint failed: C(5)\n";
        return false;
    }
    if (c_mod_iterative(10000ULL, 100000000ULL) != 37652224ULL) {
        std::cerr << "Checkpoint failed: C(10000) mod 1e8\n";
        return false;
    }
    if (c_mod_iterative(10000ULL, pow_u64(13ULL, 8)) != 617720485ULL) {
        std::cerr << "Checkpoint failed: C(10000) mod 13^8\n";
        return false;
    }

    // Internal consistency: closed-form modulo 13^8 matches iterative value for n=10000.
    if (c_mod_13_power_from_n_residue(10000ULL, false, 8) != 617720485ULL) {
        std::cerr << "Checkpoint failed: closed-form/iterative mismatch\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::cout << solve() << '\n';
    return 0;
}
