#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int kHexDigits = 12;
constexpr int kModBits = 4 * kHexDigits;  // 16^12 = 2^48
constexpr u64 kMod = 1ULL << kModBits;
constexpr u64 kMask = kMod - 1;

// We split odd numbers by residue classes modulo 2^16.
// Because 3 * 16 = 48, the class product expansion truncates after degree 2.
constexpr int kStepBits = 16;
constexpr u64 kStep = 1ULL << kStepBits;
constexpr u64 kMask32 = (1ULL << (kModBits - kStepBits)) - 1;      // 2^32 - 1
constexpr u64 kMask16 = (1ULL << (kModBits - 2 * kStepBits)) - 1;  // 2^16 - 1
static_assert(3 * kStepBits >= kModBits);

inline u64 mul_mod(u64 a, u64 b) { return static_cast<u64>((static_cast<u128>(a) * b) & kMask); }

u64 pow_mod(u64 base, u64 exp) {
    u64 result = 1;
    base &= kMask;
    while (exp > 0) {
        if (exp & 1ULL) result = mul_mod(result, base);
        base = mul_mod(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 inverse_odd_mod_2_48(u64 odd) {
    // Newton iteration x <- x * (2 - a*x) in Z/2^kZ
    u64 x = 1;
    for (int i = 0; i < 6; ++i) {
        const u64 ax = mul_mod(odd, x);
        const u64 two_minus_ax = (2 + kMod - ax) & kMask;
        x = mul_mod(x, two_minus_ax);
    }
    return x;
}

std::vector<u64> build_inverse_table() {
    std::vector<u64> inv(kStep / 2);
    for (u64 r = 1; r < kStep; r += 2) {
        inv[r >> 1] = inverse_odd_mod_2_48(r);
    }
    return inv;
}

const std::vector<u64>& odd_inverse_table() {
    static const std::vector<u64> table = build_inverse_table();
    return table;
}

u64 e1_mod_2_32(u64 c) {
    // e1 = sum_{j=0}^{c-1} j = c(c-1)/2  (mod 2^32)
    if (c < 2) return 0;
    u64 a = c;
    u64 b = c - 1;
    if ((a & 1ULL) == 0) {
        a >>= 1ULL;
    } else {
        b >>= 1ULL;
    }
    return static_cast<u64>((static_cast<u128>(a & kMask32) * (b & kMask32)) & kMask32);
}

u64 e2_mod_2_16(u64 c) {
    // e2 = sum_{0<=i<j<c} i*j = c(c-1)(c-2)(3c-1)/24  (mod 2^16)
    if (c < 3) return 0;

    std::array<u64, 4> f = {c, c - 1, c - 2, 3 * c - 1};

    int twos_to_remove = 3;  // divide by 8
    for (u64& v : f) {
        while (twos_to_remove > 0 && (v & 1ULL) == 0) {
            v >>= 1ULL;
            --twos_to_remove;
        }
    }
    assert(twos_to_remove == 0);

    bool divided_by_three = false;
    for (int i = 0; i < 3; ++i) {
        if (f[i] % 3 == 0) {
            f[i] /= 3;
            divided_by_three = true;
            break;
        }
    }
    assert(divided_by_three);

    u64 res = 1;
    for (u64 v : f) {
        res = static_cast<u64>((static_cast<u128>(res) * (v & kMask16)) & kMask16);
    }
    return res;
}

u64 class_product(u64 odd_residue, u64 count, const std::vector<u64>& inv) {
    if (count == 0) return 1;

    const u64 r_pow = pow_mod(odd_residue, count);
    const u64 r_inv = inv[odd_residue >> 1];
    const u64 a = mul_mod(kStep, r_inv);  // 2^16 / r in Z/2^48Z

    // Product_j (r + j*2^16) = r^c * Product_j (1 + j*a)
    // and (j*a)^3 vanishes mod 2^48, so keep up to degree 2.
    u64 poly = 1;
    poly = (poly + mul_mod(a, e1_mod_2_32(count))) & kMask;
    poly = (poly + mul_mod(mul_mod(a, a), e2_mod_2_16(count))) & kMask;

    return mul_mod(r_pow, poly);
}

u64 odd_prefix(u64 n, const std::vector<u64>& inv) {
    // Product of odd integers <= n modulo 2^48.
    if (n == 0) return 1;

    const u64 q = n >> kStepBits;
    const u64 rem = n & (kStep - 1);

    if (q == 0) {
        u64 p = 1;
        for (u64 r = 1; r <= rem; r += 2) {
            p = mul_mod(p, r);
        }
        return p;
    }

    u64 p = 1;
    const u64 q_step = (q * kStep) & kMask;
    for (u64 r = 1; r < kStep; r += 2) {
        p = mul_mod(p, class_product(r, q, inv));
        if (r <= rem) {
            // class(r, q+1) = class(r, q) * (r + q*2^16)
            p = mul_mod(p, (q_step + r) & kMask);
        }
    }
    return p;
}

u64 odd_factorial_part(u64 n) {
    // odd(n!) = Product_{i>=0} Product_{odd <= floor(n/2^i)} odd
    const auto& inv = odd_inverse_table();

    std::vector<u64> levels;
    for (u64 x = n; x > 0; x >>= 1ULL) {
        levels.push_back(x);
    }

    std::vector<u64> pref(levels.size(), 1);
    unsigned workers = std::thread::hardware_concurrency();
    if (workers == 0) workers = 4;
    workers = std::min<unsigned>(workers, static_cast<unsigned>(levels.size()));

    if (workers <= 1 || levels.size() < 4) {
        for (size_t i = 0; i < levels.size(); ++i) {
            pref[i] = odd_prefix(levels[i], inv);
        }
    } else {
        std::atomic<size_t> next{0};
        std::vector<std::thread> threads;
        threads.reserve(workers);

        for (unsigned t = 0; t < workers; ++t) {
            threads.emplace_back([&]() {
                while (true) {
                    const size_t idx = next.fetch_add(1);
                    if (idx >= levels.size()) break;
                    pref[idx] = odd_prefix(levels[idx], inv);
                }
            });
        }
        for (auto& th : threads) th.join();
    }

    u64 res = 1;
    for (u64 v : pref) {
        res = mul_mod(res, v);
    }
    return res;
}

u64 v2_factorial(u64 n) {
    u64 e = 0;
    while (n > 0) {
        n >>= 1ULL;
        e += n;
    }
    return e;
}

u64 f_value(u64 n) {
    const u64 odd_part = odd_factorial_part(n);
    const u64 rem_twos = v2_factorial(n) & 3ULL;
    return mul_mod(odd_part, (1ULL << rem_twos));
}

u64 brute_f_value(u64 n) {
    u64 odd_part = 1;
    u64 twos = 0;
    for (u64 i = 1; i <= n; ++i) {
        u64 x = i;
        while ((x & 1ULL) == 0) {
            x >>= 1ULL;
            ++twos;
        }
        odd_part = mul_mod(odd_part, x & kMask);
    }
    return mul_mod(odd_part, (1ULL << (twos & 3ULL)));
}

u64 factorial_u64(unsigned n) {
    u64 f = 1;
    for (unsigned i = 2; i <= n; ++i) f *= i;
    return f;
}

std::string to_hex12(u64 x) {
    std::ostringstream out;
    out << std::uppercase << std::hex << std::setw(kHexDigits) << std::setfill('0') << x;
    return out.str();
}

bool validation_checkpoints() {
    bool ok = true;

    // Provided in the statement: 20! = ...21C3677C82B40000 => f(20) = 21C3677C82B4.
    const u64 expected_20 = 0x21C3677C82B4ULL;
    const u64 got_20 = f_value(20);
    std::cout << "Check f(20):      " << to_hex12(got_20) << " (expected 21C3677C82B4)\n";
    if (got_20 != expected_20) ok = false;

    // Cross-check against brute force on smaller n values.
    const std::vector<u64> tests = {1, 2, 3, 7, 10, 25, 100, 1000, 100000};
    for (u64 n : tests) {
        const u64 fast = f_value(n);
        const u64 brute = brute_f_value(n);
        std::cout << "Cross-check n=" << n << ": " << to_hex12(fast);
        if (fast == brute) {
            std::cout << " (ok)\n";
        } else {
            std::cout << " (mismatch, brute " << to_hex12(brute) << ")\n";
            ok = false;
        }
    }

    return ok;
}

}  // namespace

int main() {
    if (!validation_checkpoints()) {
        std::cerr << "Validation failed.\n";
        return 1;
    }

    const u64 n = factorial_u64(20);
    const u64 ans = f_value(n);

    std::cout << "N = 20! = " << n << '\n';
    std::cout << "f(20!) = " << to_hex12(ans) << '\n';
    return 0;
}

