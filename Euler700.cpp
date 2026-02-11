#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kA = 1'504'170'715'041'707ULL;
constexpr u64 kM = 4'503'599'627'370'517ULL;
constexpr u64 kForwardStop = 20'000'000ULL;

u64 mod_mul(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) %
                            static_cast<u128>(mod));
}

u64 mod_inverse(u64 a, u64 mod) {
    std::int64_t t = 0;
    std::int64_t new_t = 1;
    std::int64_t r = static_cast<std::int64_t>(mod);
    std::int64_t new_r = static_cast<std::int64_t>(a);

    while (new_r != 0) {
        const std::int64_t q = r / new_r;

        const std::int64_t t_next = t - q * new_t;
        t = new_t;
        new_t = t_next;

        const std::int64_t r_next = r - q * new_r;
        r = new_r;
        new_r = r_next;
    }

    assert(r == 1);
    if (t < 0) {
        t += static_cast<std::int64_t>(mod);
    }
    return static_cast<u64>(t);
}

u64 sum_all_eulercoins() {
    u64 sum = 0;

    u64 term = kA;
    u64 best_coin = kA;
    sum += best_coin;

    while (best_coin > kForwardStop) {
        term += kA;
        if (term >= kM) {
            term -= kM;
        }
        if (term < best_coin) {
            best_coin = term;
            sum += best_coin;
        }
    }

    const u64 inv = mod_inverse(kA, kM);
    assert(mod_mul(kA, inv, kM) == 1ULL);

    u64 best_index = kM + 1ULL;
    for (u64 value = 1ULL; value < best_coin; ++value) {
        const u64 idx = mod_mul(inv, value, kM);
        if (idx < best_index) {
            best_index = idx;
            sum += value;
        }
    }

    return sum;
}

}  // namespace

int main() {
    const u64 second_eulercoin = mod_mul(kA, 3ULL, kM);
    assert(kA + second_eulercoin == 1'513'083'232'796'311ULL);

    std::cout << sum_all_eulercoins() << '\n';
    return 0;
}

