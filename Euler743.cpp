#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'000'000'007ULL;

inline u64 mod_mul(const u64 a, const u64 b) {
    return (a * b) % kMod;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        exp >>= 1ULL;
    }
    return result;
}

u64 A(const u64 k, const u64 n) {
    assert(n % k == 0ULL);

    const u64 q = n / k;
    const u64 b = mod_pow(2ULL, q);

    u64 term = mod_pow(b, k);
    u64 answer = term;

    const u64 inv_b = mod_pow(b, kMod - 2ULL);
    const u64 inv_b2 = mod_mul(inv_b, inv_b);

    const u64 half = k / 2ULL;
    std::vector<u32> inverses(static_cast<std::size_t>(half + 1ULL), 0U);
    if (half >= 1ULL) {
        inverses[1] = 1U;
    }
    for (u64 i = 2ULL; i <= half; ++i) {
        inverses[static_cast<std::size_t>(i)] = static_cast<u32>(
            kMod - mod_mul(kMod / i, inverses[static_cast<std::size_t>(kMod % i)]));
    }

    u64 num1 = k % kMod;
    u64 num2 = (k + kMod - 1ULL) % kMod;
    for (u64 x = 1ULL; x <= half; ++x) {
        const u64 inv_x = inverses[static_cast<std::size_t>(x)];

        term = mod_mul(term, num1);
        term = mod_mul(term, num2);
        term = mod_mul(term, inv_x);
        term = mod_mul(term, inv_x);
        term = mod_mul(term, inv_b2);

        answer += term;
        if (answer >= kMod) {
            answer -= kMod;
        }

        num1 += kMod - 2ULL;
        if (num1 >= kMod) {
            num1 -= kMod;
        }
        num2 += kMod - 2ULL;
        if (num2 >= kMod) {
            num2 -= kMod;
        }
    }

    return answer;
}

}  // namespace

int main() {
    assert(A(3ULL, 9ULL) == 560ULL);
    assert(A(4ULL, 20ULL) == 1'060'870ULL);

    std::cout << A(100'000'000ULL, 10'000'000'000'000'000ULL) << '\n';
    return 0;
}
