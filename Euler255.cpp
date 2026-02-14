#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

static inline u64 ceil_div(u64 a, u64 b) {
    return (a + b - 1) / b;
}

static u128 total_iterations(u64 lo, u64 hi, u64 x) {
    if (lo > hi) return 0;

    u128 total = static_cast<u128>(hi - lo + 1);  // Current step always counts.

    const u64 cmin = ceil_div(lo, x);
    const u64 cmax = ceil_div(hi, x);
    const u64 ymin = (x + cmin) / 2;
    const u64 ymax = (x + cmax) / 2;

    for (u64 y = ymin; y <= ymax; ++y) {
        i64 cl = static_cast<i64>(2 * y) - static_cast<i64>(x);
        i64 ch = cl + 1;
        i64 from = std::max<i64>(cl, static_cast<i64>(cmin));
        i64 to = std::min<i64>(ch, static_cast<i64>(cmax));
        if (from > to) continue;

        u64 nl = static_cast<u64>(from - 1) * x + 1;
        u64 nr = static_cast<u64>(to) * x;
        if (nl < lo) nl = lo;
        if (nr > hi) nr = hi;

        if (y == x) continue;
        total += total_iterations(nl, nr, y);
    }
    return total;
}

static long double average_iterations_for_digits(int d) {
    u64 lo = 1;
    for (int i = 1; i < d; ++i) lo *= 10;
    u64 hi = lo * 10 - 1;

    u64 p10 = 1;
    int exp = (d % 2 == 1) ? ((d - 1) / 2) : ((d - 2) / 2);
    for (int i = 0; i < exp; ++i) p10 *= 10;
    u64 x0 = (d % 2 == 1) ? (2 * p10) : (7 * p10);
    u128 tot = total_iterations(lo, hi, x0);
    u64 cnt = hi - lo + 1;
    return static_cast<long double>(tot) / static_cast<long double>(cnt);
}

int main() {
    long double check = average_iterations_for_digits(5);
    assert(check > 3.2102888888L && check < 3.2102888890L);

    long double ans = average_iterations_for_digits(14);
    std::cout << std::fixed << std::setprecision(10) << static_cast<double>(ans) << '\n';
    return 0;
}
