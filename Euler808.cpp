#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static inline u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(mod));
}

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    a %= mod;
    while (e > 0) {
        if (e & 1ULL) {
            r = mul_mod(r, a, mod);
        }
        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool is_prime(u64 n) {
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
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    static constexpr u64 WITNESSES[] = {2ULL, 325ULL, 9'375ULL, 28'178ULL, 450'775ULL, 9'780'504ULL, 1'795'265'022ULL};
    for (u64 a : WITNESSES) {
        if (a % n == 0) {
            continue;
        }
        u64 x = pow_mod(a, d, n);
        if (x == 1 || x == n - 1) {
            continue;
        }
        bool comp = true;
        for (int r = 1; r < s; ++r) {
            x = mul_mod(x, x, n);
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

static u64 reverse_digits(u64 x) {
    u64 r = 0;
    while (x > 0) {
        r = r * 10 + (x % 10);
        x /= 10;
    }
    return r;
}

static bool is_palindrome(u64 x) {
    return x == reverse_digits(x);
}

static u64 isqrt_u64(u64 n) {
    u64 lo = 0;
    u64 hi = 1;
    while (hi <= n / hi) {
        hi <<= 1ULL;
    }
    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;
        if (mid <= n / mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

static std::vector<u64> first_reversible_prime_squares(int need) {
    std::vector<u64> vals;
    vals.reserve(need);

    u64 p = 2;
    while (static_cast<int>(vals.size()) < need) {
        if (is_prime(p)) {
            const u64 sq = p * p;
            if (!is_palindrome(sq)) {
                const u64 rev = reverse_digits(sq);
                const u64 r = isqrt_u64(rev);
                if (r * r == rev && is_prime(r)) {
                    vals.push_back(sq);
                }
            }
        }
        if (p == 2) {
            p = 3;
        } else {
            p += 2;
        }
    }

    return vals;
}

int main() {
    const auto vals = first_reversible_prime_squares(50);

    assert(vals[0] == 169ULL);
    assert(vals[1] == 961ULL);

    u64 sum = 0;
    for (u64 v : vals) {
        sum += v;
    }

    std::cout << sum << '\n';
    return 0;
}
