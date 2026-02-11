#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

bool power_leq(const u64 base, const int exp, const u64 limit) {
    u128 p = 1;
    for (int i = 0; i < exp; ++i) {
        p *= static_cast<u128>(base);
        if (p > static_cast<u128>(limit)) {
            return false;
        }
    }
    return true;
}

u64 int_root(const u64 n, const int k) {
    long double guess = std::powl(static_cast<long double>(n), 1.0L / static_cast<long double>(k));
    u64 r = static_cast<u64>(guess);
    if (r < 1ULL) {
        r = 1ULL;
    }
    while (power_leq(r + 1ULL, k, n)) {
        ++r;
    }
    while (!power_leq(r, k, n)) {
        --r;
    }
    return r;
}

u64 count_nondecreasing(const u64 limit, const u64 min_factor, const int parts) {
    if (parts == 0) {
        return 1ULL;
    }
    if (parts == 1) {
        if (limit < min_factor) {
            return 0ULL;
        }
        return limit - min_factor + 1ULL;
    }

    const u64 r = int_root(limit, parts);
    if (r < min_factor) {
        return 0ULL;
    }

    u64 total = 0ULL;
    for (u64 x = min_factor; x <= r; ++x) {
        total += count_nondecreasing(limit / x, x, parts - 1);
    }
    return total;
}

u64 D(const u64 N, const u64 K) {
    int max_parts = 0;
    for (u64 p = 1ULL; p <= N; p <<= 1U) {
        ++max_parts;
    }
    --max_parts;

    std::vector<u64> counts(static_cast<std::size_t>(max_parts + 1), 0ULL);
    counts[0] = 1ULL;
    for (int m = 1; m <= max_parts; ++m) {
        counts[static_cast<std::size_t>(m)] = count_nondecreasing(N, 2ULL, m);
    }

    u64 answer = K % kMod;
    for (int m = 1; m <= max_parts; ++m) {
        if (static_cast<u64>(m) > K) {
            break;
        }
        const u64 ways = counts[static_cast<std::size_t>(m)] % kMod;
        const u64 multiplicity = (K - static_cast<u64>(m) + 1ULL) % kMod;
        answer = (answer + static_cast<u64>((static_cast<u128>(ways) * multiplicity) % kMod)) % kMod;
    }

    return answer;
}

}  // namespace

int main() {
    assert(D(10ULL, 10ULL) == 153ULL);
    assert(D(100ULL, 100ULL) == 35'384ULL);

    std::cout << D(10'000'000'000ULL, 10'000'000'000ULL) << '\n';
    return 0;
}
