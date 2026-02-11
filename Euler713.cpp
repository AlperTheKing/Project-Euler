#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using i64 = std::int64_t;

i64 t_value(const i64 n, const i64 m) {
    const i64 r = m - 1;
    const i64 q = n / r;
    const i64 s = n - r * q;
    return r * q * (q - 1) / 2 + s * q;
}

i64 l_value(const i64 n) {
    i64 total = 0;
    for (i64 r = 1; r <= n - 1; ++r) {
        const i64 q = n / r;
        const i64 s = n - r * q;
        total += r * q * (q - 1) / 2 + s * q;
    }
    return total;
}

}  // namespace

int main() {
    assert(t_value(3, 2) == 3);
    assert(t_value(8, 4) == 7);
    assert(l_value(1'000) == 3'281'346);

    std::cout << l_value(10'000'000) << '\n';
    return 0;
}
