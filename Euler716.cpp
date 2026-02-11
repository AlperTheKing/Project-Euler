#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using i64 = std::int64_t;

constexpr i64 kMod = 1'000'000'007LL;

i64 mod_norm(i64 x) {
    x %= kMod;
    if (x < 0) {
        x += kMod;
    }
    return x;
}

i64 mod_mul(const i64 a, const i64 b) {
    return static_cast<i64>((__int128)a * b % kMod);
}

i64 mod_pow2(i64 e) {
    i64 base = 2;
    i64 result = 1;
    while (e > 0) {
        if (e & 1LL) {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        e >>= 1LL;
    }
    return result;
}

i64 C_mod(const i64 h, const i64 w) {
    const i64 hw = mod_mul(h % kMod, w % kMod);

    const i64 p2h = mod_pow2(h);
    const i64 p2w = mod_pow2(w);
    const i64 p2hw = mod_pow2(h + w);

    const i64 t0 = mod_mul(9, p2hw);

    const i64 c1 = mod_norm(2 * hw - mod_mul(8, w % kMod) - 10);
    const i64 t1 = mod_mul(c1, p2h);

    const i64 c2 = mod_norm(2 * hw - mod_mul(8, h % kMod) - 10);
    const i64 t2 = mod_mul(c2, p2w);

    const i64 t3 = mod_norm(2 * hw + mod_mul(10, h % kMod) + mod_mul(10, w % kMod) + 10);

    return mod_norm(t0 + t1 + t2 + t3);
}

}  // namespace

int main() {
    assert(C_mod(3, 3) == 408);
    assert(C_mod(3, 6) == 4696);
    assert(C_mod(10, 20) == 988'971'143);

    std::cout << C_mod(10'000, 20'000) << '\n';
    return 0;
}
