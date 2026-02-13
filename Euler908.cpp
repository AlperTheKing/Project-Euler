#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 1'111'211'113LL;

i64 mod_pow(i64 a, i64 e) {
    i64 r = 1 % kMod;
    i64 x = a % kMod;
    i64 p = e;
    while (p > 0) {
        if (p & 1) {
            r = static_cast<i64>((__int128)r * x % kMod);
        }
        x = static_cast<i64>((__int128)x * x % kMod);
        p >>= 1;
    }
    return r;
}

std::vector<int> primes_up_to(int n) {
    std::vector<bool> is_prime(n + 1, true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; (i64)p * p <= n; ++p) {
        if (!is_prime[p]) {
            continue;
        }
        for (int q = p * p; q <= n; q += p) {
            is_prime[q] = false;
        }
    }
    std::vector<int> primes;
    for (int x = 2; x <= n; ++x) {
        if (is_prime[x]) {
            primes.push_back(x);
        }
    }
    return primes;
}

u64 triangular_residue_count_bruteforce(int s) {
    std::vector<char> seen(s, 0);
    int cnt = 0;
    for (int n = 1; n <= 2 * s; ++n) {
        const int r = (int)((1LL * n * (n + 1) / 2) % s);
        if (!seen[r]) {
            seen[r] = 1;
            ++cnt;
        }
    }
    return (u64)cnt;
}

u64 u_prime_power(int p, int e) {
    if (p == 2) {
        return u64{1} << e;
    }
    u64 u = (u64)(p + 1) / 2;
    for (int k = 2; k <= e; ++k) {
        const u64 d = (k % 2 == 0) ? (u64)(p - 1) : (u64)(p - 1) / 2;
        u = (u64)p * u - d;
    }
    return u;
}

u64 u_from_factorization(int s, const std::vector<int>& primes) {
    int x = s;
    u64 u = 1;
    for (int p : primes) {
        if ((i64)p * p > x) {
            break;
        }
        if (x % p != 0) {
            continue;
        }
        int e = 0;
        while (x % p == 0) {
            x /= p;
            ++e;
        }
        u *= u_prime_power(p, e);
    }
    if (x > 1) {
        u *= u_prime_power(x, 1);
    }
    return u;
}

std::vector<std::vector<std::pair<u64, int>>> prime_options(int n_limit) {
    const std::vector<int> primes = primes_up_to(2 * n_limit);
    std::vector<std::vector<std::pair<u64, int>>> options(primes.size());

    for (std::size_t i = 0; i < primes.size(); ++i) {
        const int p = primes[i];
        auto& v = options[i];
        if (p == 2) {
            u64 pe = 1;
            for (;;) {
                pe <<= 1;
                if ((int)pe > n_limit) {
                    break;
                }
                v.push_back({pe, (int)pe});
            }
            continue;
        }

        u64 pe = 1;
        u64 u = 0;
        int e = 0;
        for (;;) {
            ++e;
            pe *= (u64)p;
            if (e == 1) {
                u = (u64)(p + 1) / 2;
            } else {
                const u64 d = (e % 2 == 0) ? (u64)(p - 1) : (u64)(p - 1) / 2;
                u = (u64)p * u - d;
            }
            if ((int)u > n_limit) {
                break;
            }
            v.push_back({pe, (int)u});
        }
    }

    return options;
}

std::vector<std::pair<u64, int>> enumerate_S_and_u(int n_limit) {
    const auto options = prime_options(n_limit);

    std::vector<std::pair<u64, int>> out;
    out.reserve((std::size_t)n_limit * 5);

    std::function<void(std::size_t, u64, int)> dfs = [&](std::size_t idx, u64 s, int u) {
        if (u > n_limit) {
            return;
        }
        if (idx == options.size()) {
            out.push_back({s, u});
            return;
        }

        dfs(idx + 1, s, u);
        for (const auto& [pe, up] : options[idx]) {
            const i64 nu = 1LL * u * up;
            if (nu > n_limit) {
                break;
            }
            dfs(idx + 1, s * pe, (int)nu);
        }
    };

    dfs(0, 1, 1);
    return out;
}

i64 count_clock_sequences(int n_limit) {
    const auto su = enumerate_S_and_u(n_limit);

    std::vector<i64> inv(n_limit + 2, 0);
    inv[1] = 1;
    for (int i = 2; i <= n_limit + 1; ++i) {
        inv[i] = (kMod - (kMod / i) * inv[kMod % i] % kMod) % kMod;
    }

    std::vector<i64> rep(n_limit + 1, 0);

    for (const auto& [s, u] : su) {
        if (u > n_limit) {
            continue;
        }
        const u64 free = s - (u64)u;
        const int max_k = (int)std::min<u64>(free, (u64)(n_limit - u));

        i64 comb = 1;
        int p = u;
        for (int k = 0; k <= max_k; ++k) {
            rep[p] += comb;
            if (rep[p] >= kMod) {
                rep[p] -= kMod;
            }
            if (k == max_k) {
                break;
            }
            comb = static_cast<i64>((__int128)comb * (i64)(free - (u64)k) % kMod);
            comb = static_cast<i64>((__int128)comb * inv[k + 1] % kMod);
            ++p;
        }
    }

    std::vector<i64> exact = rep;
    for (int d = 1; d <= n_limit; ++d) {
        if (exact[d] == 0) {
            continue;
        }
        const i64 sub = exact[d];
        for (int m = d + d; m <= n_limit; m += d) {
            exact[m] -= sub;
            if (exact[m] < 0) {
                exact[m] += kMod;
            }
        }
    }

    i64 total = 0;
    for (int p = 1; p <= n_limit; ++p) {
        total += exact[p];
        if (total >= kMod) {
            total -= kMod;
        }
    }
    return total;
}

void validate() {
    const auto primes = primes_up_to(500);
    for (int s = 1; s <= 200; ++s) {
        assert(u_from_factorization(s, primes) == triangular_residue_count_bruteforce(s));
    }

    assert(count_clock_sequences(3) == 3);
    assert(count_clock_sequences(4) == 7);
    assert(count_clock_sequences(10) == 561);
}

}  // namespace

int main() {
    validate();
    std::cout << count_clock_sequences(10'000) << '\n';
    return 0;
}
