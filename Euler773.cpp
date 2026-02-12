#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using i64 = long long;

constexpr i64 MOD = 1'000'000'007LL;

bool is_prime(int n) {
    if (n < 2) {
        return false;
    }
    if (n % 2 == 0) {
        return n == 2;
    }
    for (int d = 3; static_cast<i64>(d) * d <= n; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }
    return true;
}

std::vector<int> first_primes_ending7(int k) {
    std::vector<int> out;
    for (int x = 7; static_cast<int>(out.size()) < k; x += 10) {
        if (is_prime(x)) {
            out.push_back(x);
        }
    }
    return out;
}

i64 mod_pow(i64 a, i64 e) {
    i64 r = 1;
    i64 cur = a % MOD;
    while (e > 0) {
        if (e & 1LL) {
            r = static_cast<i64>((static_cast<__int128>(r) * cur) % MOD);
        }
        cur = static_cast<i64>((static_cast<__int128>(cur) * cur) % MOD);
        e >>= 1LL;
    }
    return r;
}

i64 solve(int k) {
    const std::vector<int> primes = first_primes_ending7(k);

    i64 p_mod = 1;
    i64 phi_mod = 1;
    for (const int p : primes) {
        p_mod = static_cast<i64>((static_cast<__int128>(p_mod) * p) % MOD);
        phi_mod = static_cast<i64>((static_cast<__int128>(phi_mod) * (p - 1)) % MOD);
    }

    const int pattern[4] = {7, 1, 3, 9};
    i64 s = 0;
    i64 comb = 1;

    for (int t = 0; t <= k; ++t) {
        i64 term = static_cast<i64>((static_cast<__int128>(comb) * pattern[t & 3]) % MOD);
        if ((t & 1) == 0) {
            s += term;
            if (s >= MOD) {
                s -= MOD;
            }
        } else {
            s -= term;
            if (s < 0) {
                s += MOD;
            }
        }

        if (t < k) {
            const i64 num = k - t;
            const i64 inv = mod_pow(t + 1, MOD - 2);
            comb = static_cast<i64>((static_cast<__int128>(comb) * num) % MOD);
            comb = static_cast<i64>((static_cast<__int128>(comb) * inv) % MOD);
        }
    }

    i64 inner = (5LL * phi_mod + s) % MOD;
    return static_cast<i64>((static_cast<__int128>(p_mod) * inner) % MOD);
}

}  // namespace

int main() {
    assert(solve(3) == 76'101'452LL);
    std::cout << solve(97) << '\n';
    return 0;
}
