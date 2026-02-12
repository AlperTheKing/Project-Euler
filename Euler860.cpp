#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using i64 = std::int64_t;

static constexpr i64 MOD = 989'898'989LL;

static i64 F_value(int n) {
    const int span = 4 * n + 4;
    const int size = 2 * span + 1;
    const int shift = span;

    std::vector<i64> cur(size, 0), nxt(size, 0);
    cur[shift] = 1;

    for (int step = 1; step <= n; ++step) {
        std::fill(nxt.begin(), nxt.end(), 0);
        int lo = shift - 4 * step;
        int hi = shift + 4 * step;

        for (int i = lo; i <= hi; ++i) {
            i64 v = cur[i - 4];
            v += cur[i + 4];
            if (v >= MOD) v -= MOD;
            v += cur[i - 1];
            if (v >= MOD) v -= MOD;
            v += cur[i + 1];
            if (v >= MOD) v -= MOD;
            nxt[i] = v;
        }

        cur.swap(nxt);
    }

    return cur[shift];
}

int main() {
    assert(F_value(2) == 4);
    assert(F_value(10) == 63'594);
    std::cout << F_value(9'898) << '\n';
    return 0;
}
