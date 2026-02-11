#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kDiv = 1'000'000LL;

i64 mod_fix(const i64 x) {
    i64 y = x % kDiv;
    if (y < 0) {
        y += kDiv;
    }
    return y;
}

i64 t_exact_small(const int n) {
    std::array<i64, 64> a{};
    a[0] = 1;
    a[1] = 1;
    a[2] = 1;
    a[3] = 2;
    a[4] = 2;
    a[5] = 3;
    a[6] = 4;

    for (int i = 7; i <= n; ++i) {
        a[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i - 1)] +
                                         2 * a[static_cast<std::size_t>(i - 2)] -
                                         2 * a[static_cast<std::size_t>(i - 3)] -
                                         a[static_cast<std::size_t>(i - 4)] +
                                         a[static_cast<std::size_t>(i - 5)] +
                                         a[static_cast<std::size_t>(i - 6)] -
                                         a[static_cast<std::size_t>(i - 7)];
    }

    const i64 p = 1LL << (n / 2);
    return p - a[static_cast<std::size_t>(n)];
}

int solve() {
    std::array<i64, 7> a_mod{1, 1, 1, 2, 2, 3, 4};

    i64 pow2 = 1;  // 2^{floor(n/2)} mod 1e6 at current n
    for (int n = 1;; ++n) {
        if ((n & 1) == 0) {
            pow2 = (pow2 * 2) % kDiv;
        }

        if (n >= 7) {
            const i64 next = a_mod[static_cast<std::size_t>((n - 1) % 7)] +
                             2 * a_mod[static_cast<std::size_t>((n - 2) % 7)] -
                             2 * a_mod[static_cast<std::size_t>((n - 3) % 7)] -
                             a_mod[static_cast<std::size_t>((n - 4) % 7)] +
                             a_mod[static_cast<std::size_t>((n - 5) % 7)] +
                             a_mod[static_cast<std::size_t>((n - 6) % 7)] -
                             a_mod[static_cast<std::size_t>((n - 7) % 7)];
            a_mod[static_cast<std::size_t>(n % 7)] = mod_fix(next);
        }

        const i64 t_mod = mod_fix(pow2 - a_mod[static_cast<std::size_t>(n % 7)]);
        if (n > 42 && t_mod == 0) {
            return n;
        }
    }
}

}  // namespace

int main() {
    assert(t_exact_small(20) == 824);
    assert(t_exact_small(42) == 1'999'923);

    std::cout << solve() << '\n';
    return 0;
}

