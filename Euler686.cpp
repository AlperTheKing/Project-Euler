#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

using i64 = std::int64_t;

i64 p(int L, int n) {
    const int d = static_cast<int>(std::to_string(L).size());
    const long double lo = log10l(static_cast<long double>(L)) - static_cast<long double>(d - 1);
    const long double hi = log10l(static_cast<long double>(L + 1)) - static_cast<long double>(d - 1);
    const long double alpha = log10l(2.0L);

    long double x = 0.0L;
    int cnt = 0;
    i64 j = 0;

    while (cnt < n) {
        ++j;
        x += alpha;
        if (x >= 1.0L) {
            x -= 1.0L;
        }
        if (x >= lo && x < hi) {
            ++cnt;
        }
    }

    return j;
}

}  // namespace

int main() {
    assert(p(12, 1) == 7);
    assert(p(12, 2) == 80);
    assert(p(123, 45) == 12'710);

    std::cout << p(123, 678'910) << "\n";
    return 0;
}
