#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr int kHalfBits = 12;
constexpr int kKeyShift = 13;
constexpr int kValueBits = 24;
constexpr u64 kValueMask = (1ULL << kValueBits) - 1ULL;

inline int poly_degree(u32 p) {
    return p ? (31 - __builtin_clz(p)) : -1;
}

u32 poly_mod(u32 a, u32 b) {
    assert(b != 0);
    const int db = poly_degree(b);
    while (a != 0) {
        const int da = poly_degree(a);
        if (da < db) {
            break;
        }
        a ^= (b << (da - db));
    }
    return a;
}

u32 poly_gcd(u32 a, u32 b) {
    while (b != 0) {
        u32 r = poly_mod(a, b);
        a = b;
        b = r;
    }
    return a;
}

u32 poly_div_exact(u32 a, u32 b) {
    assert(b != 0);
    u32 q = 0;
    const int db = poly_degree(b);
    while (a != 0) {
        const int da = poly_degree(a);
        assert(da >= db);
        const int s = da - db;
        q ^= (1U << s);
        a ^= (b << s);
    }
    return q;
}

inline void split_even_odd_bits(u32 n, u32& even_part, u32& odd_part) {
    even_part = 0;
    odd_part = 0;
    for (int i = 0; i < kHalfBits; ++i) {
        even_part |= ((n >> (2 * i)) & 1U) << i;
        odd_part |= ((n >> (2 * i + 1)) & 1U) << i;
    }
}

inline u32 pack_key(u32 x, u32 y) {
    return (y << kKeyShift) | x;
}

inline u64 pack_record(u32 key, u32 value) {
    return (static_cast<u64>(key) << kValueBits) | static_cast<u64>(value);
}

inline u32 record_key(u64 record) {
    return static_cast<u32>(record >> kValueBits);
}

inline u32 record_value(u64 record) {
    return static_cast<u32>(record & kValueMask);
}

u64 count_ge(const std::vector<u32>& sorted_values, u32 threshold) {
    auto it = std::lower_bound(sorted_values.begin(), sorted_values.end(), threshold);
    return static_cast<u64>(sorted_values.end() - it);
}

u64 solve_fast(u32 N) {
    std::vector<u64> a_general;
    std::vector<u64> b_general;
    std::vector<u32> b_v_zero;
    std::vector<u32> b_u_zero;
    std::vector<u32> a_x_zero;
    std::vector<u32> a_y_zero;

    a_general.reserve(static_cast<size_t>(N));
    b_general.reserve(static_cast<size_t>(N));
    b_v_zero.reserve(5000);
    b_u_zero.reserve(5000);
    a_x_zero.reserve(5000);
    a_y_zero.reserve(5000);

    for (u32 b = 0; b <= N; ++b) {
        u32 U = 0;
        u32 V = 0;
        split_even_odd_bits(b, U, V);

        if (V == 0) {
            b_v_zero.push_back(b);
        }
        if (U == 0) {
            b_u_zero.push_back(b);
        }

        if (b == 0 || U == 0 || V == 0) {
            continue;
        }

        const u32 g = poly_gcd(U, V);
        const u32 x = poly_div_exact(V, g);
        const u32 y = poly_div_exact(U, g);
        b_general.push_back(pack_record(pack_key(x, y), b));
    }

    for (u32 a = 1; a <= N; ++a) {
        u32 X = 0;
        u32 Y = 0;
        split_even_odd_bits(a, X, Y);

        if (X == 0) {
            if (Y != 0) {
                a_x_zero.push_back(a);
            }
            continue;
        }
        if (Y == 0) {
            a_y_zero.push_back(a);
            continue;
        }

        const u32 uy = (Y << 1);
        const u32 d = poly_gcd(X, uy);
        const u32 x1 = poly_div_exact(X, d);
        const u32 y1 = poly_div_exact(uy, d);
        a_general.push_back(pack_record(pack_key(x1, y1), a));
    }

    std::sort(a_general.begin(), a_general.end());
    std::sort(b_general.begin(), b_general.end());

    u64 total = static_cast<u64>(N) + 1ULL;

    for (u32 a : a_x_zero) {
        total += count_ge(b_v_zero, a);
    }
    for (u32 a : a_y_zero) {
        total += count_ge(b_u_zero, a);
    }

    size_t i = 0;
    size_t j = 0;
    while (i < a_general.size()) {
        const u32 key = record_key(a_general[i]);
        size_t i_end = i;
        while (i_end < a_general.size() && record_key(a_general[i_end]) == key) {
            ++i_end;
        }

        while (j < b_general.size() && record_key(b_general[j]) < key) {
            const u32 skip_key = record_key(b_general[j]);
            while (j < b_general.size() && record_key(b_general[j]) == skip_key) {
                ++j;
            }
        }

        if (j < b_general.size() && record_key(b_general[j]) == key) {
            size_t j_end = j;
            while (j_end < b_general.size() && record_key(b_general[j_end]) == key) {
                ++j_end;
            }

            size_t ptr = j;
            for (size_t p = i; p < i_end; ++p) {
                const u32 a_value = record_value(a_general[p]);
                while (ptr < j_end && record_value(b_general[ptr]) < a_value) {
                    ++ptr;
                }
                total += static_cast<u64>(j_end - ptr);
            }
            j = j_end;
        }

        i = i_end;
    }

    return total;
}

u64 carryless_mul(u64 x, u64 y) {
    u64 r = 0;
    while (y != 0) {
        if (y & 1ULL) {
            r ^= x;
        }
        x <<= 1;
        y >>= 1;
    }
    return r;
}

bool has_solution_by_definition(u32 a, u32 b) {
    const u64 lhs = carryless_mul(a, a) ^ carryless_mul(2ULL, carryless_mul(a, b)) ^ carryless_mul(b, b);
    return (lhs & 0xAAAAAAAAAAAAAAAAULL) == 0ULL;
}

u64 solve_bruteforce(u32 N) {
    u64 cnt = 0;
    for (u32 a = 0; a <= N; ++a) {
        for (u32 b = a; b <= N; ++b) {
            if (has_solution_by_definition(a, b)) {
                ++cnt;
            }
        }
    }
    return cnt;
}

void run_validations() {
    assert(solve_fast(10) == 21ULL);
    for (u32 n = 0; n <= 40; ++n) {
        assert(solve_fast(n) == solve_bruteforce(n));
    }
    assert(solve_fast(100) == solve_bruteforce(100));
}

}  // namespace

int main() {
    run_validations();
    constexpr u32 kN = 10'000'000U;
    std::cout << solve_fast(kN) << '\n';
    return 0;
}
