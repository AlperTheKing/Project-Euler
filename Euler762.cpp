#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using i64 = long long;

constexpr i64 MOD = 1'000'000'000LL;

i64 C(std::uint32_t n) {
    static constexpr std::array<i64, 9> base = {
        1LL, 1LL, 2LL, 4LL, 9LL, 20LL, 46LL, 105LL, 243LL,
    };

    if (n < base.size()) {
        return base[n];
    }

    std::array<i64, 9> w = base;
    for (std::uint32_t k = 9; k <= n; ++k) {
        i64 next = 0;
        next += 2LL * w[8];
        next += 2LL * w[7];
        next -= 1LL * w[6];
        next -= 3LL * w[5];
        next -= 5LL * w[4];
        next += 4LL * w[3];
        next -= 2LL * w[2];
        next += 4LL * w[1];
        next %= MOD;
        if (next < 0) {
            next += MOD;
        }

        for (int i = 0; i < 8; ++i) {
            w[i] = w[i + 1];
        }
        w[8] = next;
    }
    return w[8];
}

}  // namespace

int main() {
    assert(C(2) == 2LL);
    assert(C(10) == 1301LL);
    assert(C(20) == 5'895'236LL);
    assert(C(100) == 125'923'036LL);

    std::cout << C(100'000) << '\n';
    return 0;
}
