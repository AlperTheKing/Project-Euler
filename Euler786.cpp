#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

using i64 = long long;
using u64 = unsigned long long;

namespace {

static i64 calc_lattice_count(i64 a, i64 b, i64 c) {
    if (a < b) {
        std::swap(a, b);
    }
    if (a > c) {
        return 0;
    }

    const i64 m = c / a;
    if (a == b) {
        return m * (m - 1) / 2;
    }

    const i64 k = (a - 1) / b;
    const i64 h = (c - a * m) / b;
    return calc_lattice_count(b, a - b * k, c - b * (k * m + h)) +
           k * m * (m - 1) / 2 +
           m * h;
}

struct MertensCache {
    std::unordered_map<i64, i64> memo;

    MertensCache() {
        memo.reserve(1 << 20);
        memo[0] = 0;
        memo[1] = 1;
    }

    i64 mertens(i64 n) {
        if (n <= 1) {
            return n;
        }

        auto it = memo.find(n);
        if (it != memo.end()) {
            return it->second;
        }

        i64 z = 1;
        const i64 f = static_cast<i64>(std::sqrt(static_cast<long double>(n)));

        const i64 lim = n / (f + 1);
        for (i64 i = 1; i <= lim; ++i) {
            const i64 c = n / i;
            if (c < n) {
                z -= mertens(c);
            }
        }

        for (i64 i = f; i >= 1; --i) {
            const i64 a = n / (i + 1) + 1;
            const i64 b = n / i;
            const i64 c = i;
            if (c < n) {
                z -= mertens(c) * (b - a + 1);
            }
        }

        memo.emplace(n, z);
        return z;
    }

    i64 mertens_not_mult_3(i64 n) {
        i64 z = 0;
        while (n > 0) {
            z += mertens(n);
            n /= 3;
        }
        return z;
    }
};

static u64 solve_alt(i64 N) {
    MertensCache mc;
    i64 z = 0;

    const i64 n = 3 * N + 6;
    const i64 f = static_cast<i64>(std::sqrt(static_cast<long double>(n)));

    const i64 lim = n / (f + 1);
    for (i64 i = 1; i <= lim; ++i) {
        const i64 a = i;
        const i64 b = i;
        const i64 g = n / i;
        const i64 q = calc_lattice_count(18, 10, g) - calc_lattice_count(18, 30, g);
        if (q == 0) {
            break;
        }
        z += (mc.mertens_not_mult_3(b) - mc.mertens_not_mult_3(a - 1)) * q;
    }

    for (i64 i = f; i >= 1; --i) {
        const i64 a = n / (i + 1) + 1;
        const i64 b = n / i;
        const i64 g = i;
        const i64 q = calc_lattice_count(18, 10, g) - calc_lattice_count(18, 30, g);
        if (q == 0) {
            break;
        }
        z += (mc.mertens_not_mult_3(b) - mc.mertens_not_mult_3(a - 1)) * q;
    }

    return static_cast<u64>(4 * z + 2);
}

}  // namespace

int main() {
    assert(solve_alt(10) == 6ULL);
    assert(solve_alt(100) == 478ULL);
    assert(solve_alt(1000) == 45790ULL);

    std::cout << solve_alt(1'000'000'000LL) << '\n';
    return 0;
}
