#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'000'000'007ULL;

inline u64 mul_mod(const u64 a, const u64 b) {
    return (a * b) % kMod;
}

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = mul_mod(result, base);
        }
        base = mul_mod(base, base);
        exp >>= 1ULL;
    }
    return result;
}

std::pair<u64, u64> fib_pair(const u64 n) {
    if (n == 0ULL) {
        return {0ULL, 1ULL};
    }
    const auto [a, b] = fib_pair(n >> 1ULL);

    const u64 two_b = (2ULL * b) % kMod;
    const u64 two_b_minus_a = (two_b + kMod - a) % kMod;
    const u64 c = mul_mod(a, two_b_minus_a);
    const u64 d = (mul_mod(a, a) + mul_mod(b, b)) % kMod;

    if ((n & 1ULL) == 0ULL) {
        return {c, d};
    }
    return {d, (c + d) % kMod};
}

u64 lucas(const u64 n) {
    if (n == 0ULL) {
        return 2ULL;
    }
    const auto [fn, fn1] = fib_pair(n);
    u64 value = (2ULL * fn1) % kMod;
    value = (value + kMod - fn) % kMod;
    return value;
}

void batch_invert(const std::vector<u32>& values,
                  std::vector<u32>& inverses,
                  std::vector<u32>& prefix) {
    const std::size_t m = values.size();
    inverses.resize(m);
    if (m == 0U) {
        return;
    }

    prefix.resize(m);
    prefix[0] = values[0];
    for (std::size_t i = 1; i < m; ++i) {
        prefix[i] = static_cast<u32>(mul_mod(prefix[i - 1], values[i]));
    }

    u64 suffix_inv = mod_pow(prefix[m - 1], kMod - 2ULL);
    for (std::size_t i = m; i-- > 0U;) {
        const u64 left = (i == 0U) ? 1ULL : prefix[i - 1U];
        inverses[i] = static_cast<u32>(mul_mod(suffix_inv, left));
        suffix_inv = mul_mod(suffix_inv, values[i]);
    }
}

u64 f(const u64 n) {
    const u64 r = n - 1ULL;

    u64 t = r;
    u64 coeff = 1ULL;  // T_r = 1

    u64 l_t = lucas(r);
    u64 l_t_plus_1 = lucas(r + 1ULL);

    u64 answer = 0ULL;

    constexpr u64 kBlock = 1'000'000ULL;
    std::vector<u32> denoms;
    std::vector<u32> inv_denoms;
    std::vector<u32> prefix;
    denoms.reserve(kBlock);
    inv_denoms.reserve(kBlock);
    prefix.reserve(kBlock);

    while (t >= 1ULL) {
        const u64 hi = t;
        const u64 lo = (hi > kBlock) ? (hi - kBlock + 1ULL) : 1ULL;
        const std::size_t m = static_cast<std::size_t>(hi - lo + 1ULL);

        denoms.resize(m);
        for (std::size_t i = 0; i < m; ++i) {
            const u64 cur_t = hi - static_cast<u64>(i);
            if (cur_t == 1ULL) {
                denoms[i] = 1ULL;
            } else {
                const u64 u = r - cur_t + 1ULL;
                denoms[i] = static_cast<u32>(mul_mod(cur_t % kMod, u % kMod));
            }
        }
        batch_invert(denoms, inv_denoms, prefix);

        for (std::size_t i = 0; i < m; ++i) {
            const u64 cur_t = hi - static_cast<u64>(i);

            answer += mul_mod(coeff, l_t_plus_1);
            if (answer >= kMod) {
                answer -= kMod;
            }

            if (cur_t == 1ULL) {
                break;
            }

            const u64 num_a = cur_t - 1ULL;
            const u64 num_b = 2ULL * r - cur_t;
            const u64 num = mul_mod(num_a % kMod, num_b % kMod);

            coeff = mul_mod(coeff, num);
            coeff = mul_mod(coeff, inv_denoms[i]);

            const u64 next_l_t_plus_1 = l_t;
            const u64 next_l_t = (l_t_plus_1 + kMod - l_t) % kMod;
            l_t_plus_1 = next_l_t_plus_1;
            l_t = next_l_t;
        }

        if (lo == 1ULL) {
            break;
        }
        t = lo - 1ULL;
    }

    return answer;
}

}  // namespace

int main() {
    assert(f(8ULL) == 2663ULL);
    assert(f(20ULL) == 742'296'999ULL);

    std::cout << f(100'000'000ULL) << '\n';
    return 0;
}
