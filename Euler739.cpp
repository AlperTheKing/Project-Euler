#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_pow(u64 base, u64 exp) {
    u64 result = 1ULL;
    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * base) % kMod);
        }
        base = static_cast<u64>((static_cast<u128>(base) * base) % kMod);
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
    const u64 c = static_cast<u64>((static_cast<u128>(a) * two_b_minus_a) % kMod);
    const u64 d = (static_cast<u64>((static_cast<u128>(a) * a) % kMod) +
                   static_cast<u64>((static_cast<u128>(b) * b) % kMod)) % kMod;

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

void batch_invert(const std::vector<u64>& values, std::vector<u64>& inverses) {
    const std::size_t m = values.size();
    inverses.assign(m, 0ULL);
    if (m == 0U) {
        return;
    }

    std::vector<u64> prefix(m, 0ULL);
    prefix[0] = values[0];
    for (std::size_t i = 1; i < m; ++i) {
        prefix[i] = static_cast<u64>((static_cast<u128>(prefix[i - 1]) * values[i]) % kMod);
    }

    u64 suffix_inv = mod_pow(prefix[m - 1], kMod - 2ULL);
    for (std::size_t i = m; i-- > 0U;) {
        const u64 left = (i == 0U) ? 1ULL : prefix[i - 1U];
        inverses[i] = static_cast<u64>((static_cast<u128>(suffix_inv) * left) % kMod);
        suffix_inv = static_cast<u64>((static_cast<u128>(suffix_inv) * values[i]) % kMod);
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
    std::vector<u64> denoms;
    std::vector<u64> inv_denoms;
    denoms.reserve(kBlock);

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
                denoms[i] = static_cast<u64>((static_cast<u128>(cur_t % kMod) * (u % kMod)) % kMod);
            }
        }
        batch_invert(denoms, inv_denoms);

        for (std::size_t i = 0; i < m; ++i) {
            const u64 cur_t = hi - static_cast<u64>(i);

            answer += static_cast<u64>((static_cast<u128>(coeff) * l_t_plus_1) % kMod);
            if (answer >= kMod) {
                answer -= kMod;
            }

            if (cur_t == 1ULL) {
                break;
            }

            const u64 num_a = cur_t - 1ULL;
            const u64 num_b = 2ULL * r - cur_t;
            const u64 num = static_cast<u64>((static_cast<u128>(num_a % kMod) * (num_b % kMod)) % kMod);

            coeff = static_cast<u64>((static_cast<u128>(coeff) * num) % kMod);
            coeff = static_cast<u64>((static_cast<u128>(coeff) * inv_denoms[i]) % kMod);

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
