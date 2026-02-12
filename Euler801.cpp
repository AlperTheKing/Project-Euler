#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i64 = std::int64_t;

static constexpr i64 MOD = 993'353'399LL;

static inline u64 mul_mod_u64(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

static u64 pow_mod_u64(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod_u64(r, a, mod);
        }
        a = mul_mod_u64(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool is_prime_u64(u64 n) {
    if (n < 2) {
        return false;
    }
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0) {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    static constexpr u64 WITNESSES[] = {2ULL, 325ULL, 9'375ULL, 28'178ULL, 450'775ULL, 9'780'504ULL, 1'795'265'022ULL};

    for (u64 a : WITNESSES) {
        if (a % n == 0) {
            continue;
        }
        u64 x = pow_mod_u64(a, d, n);
        if (x == 1 || x == n - 1) {
            continue;
        }
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod_u64(x, x, n);
            if (x == n - 1) {
                comp = false;
                break;
            }
        }
        if (comp) {
            return false;
        }
    }
    return true;
}

static u64 pollard_rho(u64 n, std::mt19937_64& rng) {
    if ((n & 1ULL) == 0ULL) {
        return 2;
    }

    std::uniform_int_distribution<u64> dist(2, n - 2);

    while (true) {
        const u64 c = dist(rng);
        u64 x = dist(rng);
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) {
            return (mul_mod_u64(v, v, n) + c) % n;
        };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            const u64 diff = x > y ? x - y : y - x;
            d = std::gcd(diff, n);
        }

        if (d != n) {
            return d;
        }
    }
}

static void factor_rec(u64 n, std::vector<u64>& out, std::mt19937_64& rng) {
    if (n == 1) {
        return;
    }
    if (is_prime_u64(n)) {
        out.push_back(n);
        return;
    }
    const u64 d = pollard_rho(n, rng);
    factor_rec(d, out, rng);
    factor_rec(n / d, out, rng);
}

static i64 mod_pow_i64(i64 a, i64 e) {
    i64 r = 1 % MOD;
    a %= MOD;
    if (a < 0) {
        a += MOD;
    }
    while (e > 0) {
        if (e & 1LL) {
            r = static_cast<i64>((static_cast<u128>(r) * static_cast<u128>(a)) % static_cast<u128>(MOD));
        }
        a = static_cast<i64>((static_cast<u128>(a) * static_cast<u128>(a)) % static_cast<u128>(MOD));
        e >>= 1LL;
    }
    return r;
}

static i64 f_prime_mod(u64 p, std::mt19937_64& rng) {
    const u64 n = p - 1;
    std::vector<u64> fac;
    fac.reserve(16);
    factor_rec(n, fac, rng);
    std::sort(fac.begin(), fac.end());

    std::map<u64, int> expo;
    for (u64 q : fac) {
        ++expo[q];
    }

    i64 h = 1;
    for (const auto& [q, e] : expo) {
        const i64 qm = static_cast<i64>(q % static_cast<u64>(MOD));
        const i64 q_e_minus_1 = mod_pow_i64(qm, e - 1);
        const i64 q_e = mod_pow_i64(qm, e);
        const i64 q_e_plus_1 = mod_pow_i64(qm, e + 1);

        i64 t = q_e_plus_1 + q_e - 1;
        t %= MOD;
        if (t < 0) {
            t += MOD;
        }

        const i64 term = static_cast<i64>((static_cast<u128>(q_e_minus_1) * static_cast<u128>(t)) % static_cast<u128>(MOD));
        h = static_cast<i64>((static_cast<u128>(h) * static_cast<u128>(term)) % static_cast<u128>(MOD));
    }

    const i64 nm = static_cast<i64>(n % static_cast<u64>(MOD));
    const i64 n2 = static_cast<i64>((static_cast<u128>(nm) * static_cast<u128>(nm)) % static_cast<u128>(MOD));
    const i64 nh = static_cast<i64>((static_cast<u128>(nm) * static_cast<u128>(h)) % static_cast<u128>(MOD));

    i64 f = n2 + nh;
    f %= MOD;
    if (f < 0) {
        f += MOD;
    }
    return f;
}

static std::vector<u64> primes_in_range(u64 lo, u64 hi) {
    std::vector<u64> out;
    if (hi < 2 || lo > hi) {
        return out;
    }

    if (lo <= 2 && 2 <= hi) {
        out.push_back(2);
    }

    u64 start = (lo <= 3 ? 3 : lo);
    if ((start & 1ULL) == 0ULL) {
        ++start;
    }

    for (u64 x = start; x <= hi; x += 2) {
        if (is_prime_u64(x)) {
            out.push_back(x);
        }
    }

    return out;
}

static i64 S_mod(u64 M, u64 N) {
    auto primes = primes_in_range(M, N);

    std::mt19937_64 rng(0xC0FFEEULL);
    i64 ans = 0;
    for (u64 p : primes) {
        ans += f_prime_mod(p, rng);
        ans %= MOD;
    }

    if (ans < 0) {
        ans += MOD;
    }
    return ans;
}

int main() {
    assert(S_mod(1, 100) == 7'381'000LL);
    assert(S_mod(1, 100'000) == 701'331'986LL);

    std::cout << S_mod(10'000'000'000'000'000ULL, 10'000'000'000'000'000ULL + 1'000'000ULL) << '\n';
    return 0;
}
