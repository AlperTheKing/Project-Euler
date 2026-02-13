#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

int bit_length(u64 x) {
    return 64 - __builtin_clzll(x);
}

bool losing_formula(u64 a, u64 b) {
    const u64 m = std::min(a, b);
    const u64 M = std::max(a, b);
    const int k = bit_length(m);
    const u64 p = u64{1} << k;
    return (M & (p - 1)) == (p - 1);
}

u64 L_formula(u64 n) {
    u128 unordered = 0;
    for (int k = 1;; ++k) {
        const u64 lo = u64{1} << (k - 1);
        if (lo > n) {
            break;
        }
        const u64 p = u64{1} << k;
        const u64 hi = std::min<u64>(n, p - 1);
        const u64 A = hi - lo + 1;
        const u64 C = (n + 1) / p;
        unordered += static_cast<u128>(A) * static_cast<u128>(C);
    }

    u64 diagonal = 0;
    for (int k = 1;; ++k) {
        const u64 x = (u64{1} << k) - 1;
        if (x > n) {
            break;
        }
        ++diagonal;
    }

    const u128 ordered = 2 * unordered - diagonal;
    return static_cast<u64>(ordered);
}

u64 L_bruteforce(int n) {
    std::vector<std::vector<char>> lose(n + 1, std::vector<char>(n + 1, 0));

    for (int s = 2; s <= 2 * n; ++s) {
        for (int a = 1; a <= n; ++a) {
            const int b = s - a;
            if (b < 1 || b > n) {
                continue;
            }

            const int m = std::min(a, b);
            bool winning = false;
            for (int x = 0; x <= m; ++x) {
                const int y = m - x;
                if (x >= a || y >= b) {
                    continue;
                }
                if (lose[a - x][b - y]) {
                    winning = true;
                    break;
                }
            }
            lose[a][b] = static_cast<char>(!winning);
        }
    }

    u64 cnt = 0;
    for (int a = 1; a <= n; ++a) {
        for (int b = 1; b <= n; ++b) {
            cnt += static_cast<u64>(lose[a][b] != 0);
            assert(lose[a][b] == static_cast<char>(losing_formula(a, b)));
        }
    }
    return cnt;
}

void validate() {
    assert(L_bruteforce(7) == 21);
    assert(L_bruteforce(49) == 221);
    assert(L_formula(7) == 21);
    assert(L_formula(49) == 221);
}

}  // namespace

int main() {
    validate();
    u64 n = 1;
    for (int i = 0; i < 17; ++i) {
        n *= 7;
    }
    std::cout << L_formula(n) << '\n';
    return 0;
}
