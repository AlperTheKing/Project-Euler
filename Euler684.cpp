#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;

constexpr i64 MOD = 1'000'000'007LL;

i64 mod_pow(i64 a, u64 e) {
    i64 r = 1;
    i64 x = a % MOD;
    u64 p = e;
    while (p > 0) {
        if ((p & 1ULL) != 0ULL) {
            r = static_cast<i64>((static_cast<i128>(r) * x) % MOD);
        }
        p >>= 1ULL;
        if (p > 0) {
            x = static_cast<i64>((static_cast<i128>(x) * x) % MOD);
        }
    }
    return r;
}

i64 norm(i64 v) {
    v %= MOD;
    if (v < 0) {
        v += MOD;
    }
    return v;
}

i64 S_mod(u64 k) {
    const u64 q = k / 9ULL;
    const u64 r = k % 9ULL;

    const i64 pow10 = mod_pow(10LL, q);

    i64 base = 6LL * (pow10 - 1LL);
    base = norm(base);
    base = norm(base - 9LL * static_cast<i64>(q % static_cast<u64>(MOD)));

    const i64 coeff = static_cast<i64>((r * (r + 3ULL) / 2ULL) % static_cast<u64>(MOD));
    i64 part = static_cast<i64>((static_cast<i128>(coeff) * pow10) % MOD);
    part = norm(part - static_cast<i64>(r));

    return norm(base + part);
}

u64 s_small(int n) {
    if (n == 0) {
        return 0;
    }
    int sum = n;
    std::vector<int> digits;
    while (sum > 9) {
        digits.push_back(9);
        sum -= 9;
    }
    digits.push_back(sum);
    u64 out = 0;
    for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
        out = out * 10ULL + static_cast<u64>(digits[static_cast<std::size_t>(i)]);
    }
    return out;
}

u64 S_small(int k) {
    u64 out = 0;
    for (int i = 1; i <= k; ++i) {
        out += s_small(i);
    }
    return out;
}

}  // namespace

int main() {
    assert(S_small(20) == 1'074ULL);
    assert(S_mod(20ULL) == 1'074LL);

    std::vector<u64> fib(91, 0ULL);
    fib[0] = 0ULL;
    fib[1] = 1ULL;
    for (int i = 2; i <= 90; ++i) {
        fib[static_cast<std::size_t>(i)] = fib[static_cast<std::size_t>(i - 1)] +
                                           fib[static_cast<std::size_t>(i - 2)];
    }

    i64 ans = 0;
    for (int i = 2; i <= 90; ++i) {
        ans = norm(ans + S_mod(fib[static_cast<std::size_t>(i)]));
    }

    std::cout << ans << "\n";
    return 0;
}
