#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;

constexpr i64 MOD = 1'000'000'007LL;
constexpr i64 INV2 = (MOD + 1LL) / 2LL;

i64 norm(i64 x) {
    x %= MOD;
    if (x < 0) {
        x += MOD;
    }
    return x;
}

i64 mul_mod(i64 a, i64 b) {
    return static_cast<i64>((static_cast<i128>(a) * static_cast<i128>(b)) % MOD);
}

i64 S_mod(u64 N) {
    i64 ans = 0;

    for (u64 k = 1;; ++k) {
        const u64 c = k * (k - 1ULL) / 2ULL;
        if (c >= N) {
            break;
        }

        const u64 M = N - c;
        const u64 q = M / k;
        const u64 r = M % k;

        i64 term1 = mul_mod(static_cast<i64>(k % static_cast<u64>(MOD)),
                            mul_mod(static_cast<i64>(q % static_cast<u64>(MOD)),
                                    static_cast<i64>((q - 1ULL) % static_cast<u64>(MOD))));
        term1 = mul_mod(term1, INV2);

        i64 term2 = mul_mod(static_cast<i64>(q % static_cast<u64>(MOD)),
                            static_cast<i64>((r + 1ULL) % static_cast<u64>(MOD)));

        ans = norm(ans + term1 + term2);
    }

    return ans;
}

u64 f_small(u64 n, u64 k) {
    if (k == 0ULL) {
        return 0ULL;
    }
    if (n < k * (k + 1ULL) / 2ULL) {
        return 0ULL;
    }
    return (n - k * (k - 1ULL) / 2ULL) / k;
}

u64 F_small(u64 n) {
    u64 out = 0ULL;
    for (u64 k = 1ULL;; ++k) {
        const u64 v = f_small(n, k);
        if (v == 0ULL && k * (k + 1ULL) / 2ULL > n) {
            break;
        }
        out += v;
    }
    return out;
}

u64 S_small(u64 N) {
    u64 out = 0ULL;
    for (u64 n = 1ULL; n <= N; ++n) {
        out += F_small(n);
    }
    return out;
}

}  // namespace

int main() {
    assert(S_small(100ULL) == 12'656ULL);
    assert(S_mod(100ULL) == 12'656LL);

    std::cout << S_mod(10'000'000'000'000'000ULL) << "\n";
    return 0;
}
