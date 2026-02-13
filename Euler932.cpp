#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

u64 pow10_int(int e) {
    u64 v = 1;
    for (int i = 0; i < e; ++i) {
        v *= 10;
    }
    return v;
}

std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(limit + 1, false);
    std::vector<int> primes;
    primes.reserve(limit / 10);

    for (int i = 2; i <= limit; ++i) {
        if (!composite[i]) {
            primes.push_back(i);
            if (static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(limit)) {
                for (int j = i * i; j <= limit; j += i) {
                    composite[j] = true;
                }
            }
        }
    }

    return primes;
}

std::vector<u64> factor_prime_powers(u64 x, const std::vector<int>& primes) {
    std::vector<u64> parts;

    for (int p : primes) {
        const u64 pp = static_cast<u64>(p);
        if (static_cast<u128>(pp) * pp > x) {
            break;
        }
        if (x % pp != 0) {
            continue;
        }
        u64 pk = 1;
        while (x % pp == 0) {
            x /= pp;
            pk *= pp;
        }
        parts.push_back(pk);
    }

    if (x > 1) {
        parts.push_back(x);
    }

    return parts;
}

i128 egcd(i128 a, i128 b, i128& x, i128& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    i128 x1 = 0;
    i128 y1 = 0;
    i128 g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

u64 mod_inv(u64 a, u64 mod) {
    i128 x = 0;
    i128 y = 0;
    i128 g = egcd(static_cast<i128>(a), static_cast<i128>(mod), x, y);
    assert(g == 1);
    i128 r = x % static_cast<i128>(mod);
    if (r < 0) {
        r += static_cast<i128>(mod);
    }
    return static_cast<u64>(r);
}

std::pair<u64, u64> crt_merge(u64 r1, u64 m1, u64 r2, u64 m2) {
    if (m1 == 1) {
        return {r2 % m2, m2};
    }
    const u64 inv = mod_inv(m1 % m2, m2);

    i128 diff = static_cast<i128>(r2) - static_cast<i128>(r1);
    diff %= static_cast<i128>(m2);
    if (diff < 0) {
        diff += static_cast<i128>(m2);
    }

    const u64 t = static_cast<u64>((diff * static_cast<i128>(inv)) % static_cast<i128>(m2));
    const u128 mod = static_cast<u128>(m1) * static_cast<u128>(m2);
    const u128 r = static_cast<u128>(r1) + static_cast<u128>(m1) * static_cast<u128>(t);
    return {static_cast<u64>(r % mod), static_cast<u64>(mod)};
}

u64 digit_count(u64 x) {
    u64 d = 0;
    do {
        ++d;
        x /= 10;
    } while (x > 0);
    return d;
}

bool is_2025_number(u64 n, const std::vector<u64>& pow10) {
    const int digits = static_cast<int>(digit_count(n));

    for (int d = 1; d < digits; ++d) {
        const u64 p10 = pow10[d];
        const u64 a = n / p10;
        const u64 b = n % p10;
        if (a == 0 || b == 0) {
            continue;
        }
        if (d > 1 && b < pow10[d - 1]) {
            continue;
        }
        const u64 s = a + b;
        if (static_cast<u128>(s) * s == n) {
            return true;
        }
    }

    return false;
}

u64 brute_sum_digits(int nd) {
    const u64 max_n = pow10_int(nd) - 1;
    u64 root = static_cast<u64>(std::sqrt(static_cast<long double>(max_n)));
    while (static_cast<u128>(root + 1) * (root + 1) <= max_n) {
        ++root;
    }
    while (static_cast<u128>(root) * root > max_n) {
        --root;
    }

    std::vector<u64> pow10(nd + 1, 1);
    for (int i = 1; i <= nd; ++i) {
        pow10[i] = pow10[i - 1] * 10;
    }

    u128 total = 0;
    for (u64 s = 2; s <= root; ++s) {
        const u64 n = static_cast<u64>(static_cast<u128>(s) * s);
        if (is_2025_number(n, pow10)) {
            total += n;
        }
    }

    return static_cast<u64>(total);
}

u64 solve_digits(int nd) {
    const u64 max_n = pow10_int(nd) - 1;
    u64 s_limit = static_cast<u64>(std::sqrt(static_cast<long double>(max_n)));
    while (static_cast<u128>(s_limit + 1) * (s_limit + 1) <= max_n) {
        ++s_limit;
    }
    while (static_cast<u128>(s_limit) * s_limit > max_n) {
        --s_limit;
    }

    std::vector<u64> pow10(nd + 1, 1);
    for (int i = 1; i <= nd; ++i) {
        pow10[i] = pow10[i - 1] * 10;
    }

    const int sieve_limit = static_cast<int>(std::sqrt(static_cast<long double>(pow10[nd - 1] - 1))) + 10;
    const std::vector<int> primes = sieve_primes(sieve_limit);

    std::unordered_set<u64> seen;
    seen.reserve(256);

    for (int d = 1; d < nd; ++d) {
        const u64 mod = pow10[d] - 1;
        const std::vector<u64> prime_powers = factor_prime_powers(mod, primes);

        std::vector<std::pair<u64, u64>> residues = {{0, 1}};
        for (u64 pk : prime_powers) {
            std::vector<std::pair<u64, u64>> next;
            next.reserve(residues.size() * 2);
            for (const auto& [r, m] : residues) {
                next.push_back(crt_merge(r, m, 0, pk));
                next.push_back(crt_merge(r, m, 1 % pk, pk));
            }
            residues.swap(next);
        }

        const u64 s_max = std::min(s_limit, pow10[d] - 1);
        const u64 b_low = (d == 1) ? 1 : pow10[d - 1];
        const u64 p10d = pow10[d];

        for (const auto& [r0, step] : residues) {
            u64 s = r0;
            if (s <= 1) {
                const u64 need = 2 - s;
                const u64 k = (need + step - 1) / step;
                if (k > (s_max - s) / step) {
                    continue;
                }
                s += k * step;
            }

            for (; s <= s_max; s += step) {
                const u128 ss = static_cast<u128>(s) * static_cast<u128>(s - 1);
                if (ss % mod != 0) {
                    continue;
                }

                const u64 a = static_cast<u64>(ss / mod);
                if (a == 0 || a >= s) {
                    continue;
                }
                const u64 b = s - a;

                if (b < b_low || b >= p10d) {
                    continue;
                }

                const u64 n = static_cast<u64>(static_cast<u128>(s) * s);
                if (n > max_n) {
                    continue;
                }
                if (static_cast<u128>(a) * p10d + b != n) {
                    continue;
                }

                seen.insert(n);
            }
        }
    }

    u128 total = 0;
    for (u64 v : seen) {
        total += v;
    }
    return static_cast<u64>(total);
}

void run_validations() {
    assert(brute_sum_digits(4) == 5131);
    assert(solve_digits(4) == 5131);

    for (int nd : {5, 6, 7}) {
        assert(solve_digits(nd) == brute_sum_digits(nd));
    }
}

}  // namespace

int main() {
    run_validations();
    std::cout << solve_digits(16) << '\n';
    return 0;
}
