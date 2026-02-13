#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 mul_mod_u64(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 pow_mod_u64(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    while (exp > 0) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod_u64(result, cur, mod);
        }
        cur = mul_mod_u64(cur, cur, mod);
        exp >>= 1ULL;
    }
    return result;
}

bool is_prime_u64(u64 n) {
    if (n < 2ULL) {
        return false;
    }
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0ULL) {
            return false;
        }
    }

    u64 d = n - 1ULL;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    for (u64 a : {2ULL, 325ULL, 9'375ULL, 28'178ULL, 450'775ULL, 9'780'504ULL, 1'795'265'022ULL}) {
        if (a % n == 0ULL) {
            continue;
        }
        u64 x = pow_mod_u64(a % n, d, n);
        if (x == 1ULL || x == n - 1ULL) {
            continue;
        }
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod_u64(x, x, n);
            if (x == n - 1ULL) {
                witness = false;
                break;
            }
        }
        if (witness) {
            return false;
        }
    }
    return true;
}

u64 pollard_rho(u64 n, std::mt19937_64& rng) {
    if ((n & 1ULL) == 0ULL) {
        return 2ULL;
    }
    if (n % 3ULL == 0ULL) {
        return 3ULL;
    }

    std::uniform_int_distribution<u64> dist(2ULL, n - 2ULL);
    while (true) {
        const u64 c = dist(rng);
        u64 x = dist(rng);
        u64 y = x;
        u64 d = 1ULL;

        auto f = [&](u64 v) {
            return (mul_mod_u64(v, v, n) + c) % n;
        };

        while (d == 1ULL) {
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

void factor_u64(u64 n, std::map<u64, int>& factors, std::mt19937_64& rng) {
    if (n == 1ULL) {
        return;
    }
    if (is_prime_u64(n)) {
        ++factors[n];
        return;
    }
    const u64 d = pollard_rho(n, rng);
    factor_u64(d, factors, rng);
    factor_u64(n / d, factors, rng);
}

u128 isqrt_u128(u128 x) {
    long double approx = std::sqrt(static_cast<long double>(x));
    u128 r = static_cast<u128>(approx);
    while ((r + 1) <= x / (r + 1)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u128 triangle_u128(u128 m) {
    return (m * (m + 1)) / 2;
}

std::pair<u128, u128> next_triangle_state(u128 m, std::mt19937_64& rng) {
    const u64 a = static_cast<u64>(m);
    const u64 b = static_cast<u64>(m + 1);

    std::map<u64, int> fac;
    factor_u64(a, fac, rng);
    factor_u64(b, fac, rng);

    std::vector<std::pair<u64, int>> pf(fac.begin(), fac.end());
    const u128 p = m * (m + 1);
    const u128 root = isqrt_u128(p);

    u128 best_u = 0;
    auto dfs = [&](auto&& self, std::size_t idx, u128 cur) -> void {
        if (idx == pf.size()) {
            if (cur > root) {
                return;
            }
            const u128 v = p / cur;
            if (((cur + v) & 1U) == 1U && v > cur + 1 && cur > best_u) {
                best_u = cur;
            }
            return;
        }

        const auto [prime, exp] = pf[idx];
        u128 val = 1;
        for (int e = 0; e <= exp; ++e) {
            self(self, idx + 1, cur * val);
            val *= static_cast<u128>(prime);
        }
    };

    dfs(dfs, 0, 1);
    assert(best_u > 0);

    const u128 v = p / best_u;
    const u128 jump = (v - best_u - 1) / 2;
    const u128 next_m = (best_u + v - 1) / 2;
    return {jump, next_m};
}

std::pair<u128, u128> solve_kth_triangle_hit(int k) {
    std::mt19937_64 rng(0xC0FFEEULL);

    int hit = 1;
    u128 index = 0;
    u128 m = 2;

    while (hit < k) {
        const auto [jump, next_m] = next_triangle_state(m, rng);
        index += jump;
        m = next_m;
        ++hit;
    }

    return {index, m};
}

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }
    std::string s;
    while (x > 0) {
        const int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

void run_validations() {
    const auto [idx10, m10] = solve_kth_triangle_hit(10);
    assert(idx10 == 2'964);
    assert(triangle_u128(m10) == 1'439'056);
}

}  // namespace

int main() {
    run_validations();
    const auto [idx70, _m70] = solve_kth_triangle_hit(70);
    std::cout << to_string_u128(idx70) << '\n';
    return 0;
}
