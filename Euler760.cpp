#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>

namespace {

using i64 = long long;
using i128 = __int128_t;

constexpr i64 MOD = 1'000'000'007LL;
constexpr i64 INV2 = (MOD + 1) / 2;

struct HS {
    i64 h;
    i64 s;
};

std::unordered_map<i64, HS> memo;

i64 add_mod(i64 a, i64 b) {
    i64 r = a + b;
    if (r >= MOD) {
        r -= MOD;
    }
    return r;
}

i64 mul_mod(i64 a, i64 b) {
    return static_cast<i64>((static_cast<i128>(a) * b) % MOD);
}

HS solve(i64 n) {
    if (n <= 0) {
        return HS{0, 0};
    }
    auto it = memo.find(n);
    if (it != memo.end()) {
        return it->second;
    }

    const i64 m = n / 2;
    const HS sm = solve(m);
    const HS sm1 = solve(m - 1);
    const i64 mm = m % MOD;

    i64 h = 0;
    i64 s = 0;

    if ((n & 1LL) == 0LL) {
        h = (mul_mod(2, sm.h) + mul_mod(2, sm1.h) + mm) % MOD;

        i64 poly = mul_mod(3, mul_mod(mm, (mm + 1) % MOD));
        poly = mul_mod(poly, INV2);
        s = (mul_mod(2, sm.s) + mul_mod(6, sm1.s) + poly) % MOD;
    } else {
        h = (mul_mod(4, sm.h) + mul_mod(2, (mm + 1) % MOD)) % MOD;

        i64 poly = 0;
        poly = (poly + mul_mod(3, mul_mod(mm, mm))) % MOD;
        poly = (poly + mul_mod(7, mm)) % MOD;
        poly = (poly + 4) % MOD;
        poly = mul_mod(poly, INV2);
        s = (mul_mod(6, sm.s) + mul_mod(2, sm1.s) + poly) % MOD;
    }

    const HS out{h, s};
    memo.emplace(n, out);
    return out;
}

i64 brute_G(int n) {
    i64 total = 0;
    for (int x = 0; x <= n; ++x) {
        for (int k = 0; k <= x; ++k) {
            total += 2LL * static_cast<i64>(k | (x - k));
        }
    }
    return total;
}

}  // namespace

int main() {
    assert(brute_G(10) == 754);
    assert(brute_G(100) == 583766);

    memo.reserve(256);
    const i64 g10 = mul_mod(2, solve(10).s);
    const i64 g100 = mul_mod(2, solve(100).s);
    assert(g10 == 754);
    assert(g100 == 583766);

    const i64 answer = mul_mod(2, solve(1'000'000'000'000'000'000LL).s);
    std::cout << answer << '\n';
    return 0;
}
