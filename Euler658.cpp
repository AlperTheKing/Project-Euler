#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_mul(u64 a, u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % kMod);
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1;
    base %= kMod;
    while (exp > 0) {
        if (exp & 1ULL) result = mod_mul(result, base);
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 solve_case(int k, u64 n) {
    std::vector<u32> inv(static_cast<std::size_t>(k + 2), 0);
    inv[1] = 1;
    for (int i = 2; i <= k + 1; ++i) {
        inv[static_cast<std::size_t>(i)] =
            static_cast<u32>(kMod - mod_mul(kMod / static_cast<u64>(i),
                                            inv[static_cast<std::size_t>(kMod % static_cast<u64>(i))]));
    }

    const u64 exp = n + 1ULL;
    const u64 n1 = exp % kMod;

    auto geometric_sum = [&](int r) -> u64 {
        if (r == 0) return 1;
        if (r == 1) return n1;
        const u64 p = mod_pow(static_cast<u64>(r), exp);
        const u64 num = (p + kMod - 1) % kMod;
        return mod_mul(num, inv[static_cast<std::size_t>(r - 1)]);
    };

    u64 answer = 0;

    u64 A_next = static_cast<u64>(k) % kMod;         // A_{k-1}
    u64 C_curr = static_cast<u64>(k + 1) % kMod;     // C(k+1, k)

    for (int r = k - 1; r >= 0; --r) {
        u64 A = 0;
        if (r == k - 1) {
            A = A_next;
        } else {
            const bool positive = (((k + 1 - r) & 1) == 0);
            const u64 signed_C = positive ? C_curr : (kMod - C_curr);
            A = (mod_mul(2, A_next) + signed_C + kMod - 1) % kMod;
            A_next = A;
        }

        answer = (answer + mod_mul(A, geometric_sum(r))) % kMod;

        if (r > 0) {
            C_curr = mod_mul(C_curr, static_cast<u64>(r + 1));
            C_curr = mod_mul(C_curr, inv[static_cast<std::size_t>(k - r + 1)]);
        }
    }

    return answer;
}

}  // namespace

int main() {
    assert(solve_case(4, 4) == 406);
    assert(solve_case(8, 8) == 27'902'680);
    assert(solve_case(10, 100) == 983'602'076);

    std::cout << solve_case(10'000'000, 1'000'000'000'000ULL) << "\n";
    return 0;
}
