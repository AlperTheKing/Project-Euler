#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 S(const u64 N) {
    u64 root = 1ULL;
    while ((root + 1ULL) * (root + 1ULL) <= N) {
        ++root;
    }

    std::vector<u64> j2(static_cast<std::size_t>(root + 1ULL), 0ULL);
    for (u64 i = 1ULL; i <= root; ++i) {
        j2[static_cast<std::size_t>(i)] = i * i;
    }

    std::vector<std::uint8_t> composite(static_cast<std::size_t>(root + 1ULL), 0U);
    for (u64 p = 2ULL; p <= root; ++p) {
        if (composite[static_cast<std::size_t>(p)] != 0U) {
            continue;
        }
        const u64 p2 = p * p;
        for (u64 m = p; m <= root; m += p) {
            composite[static_cast<std::size_t>(m)] = 1U;
            u64& v = j2[static_cast<std::size_t>(m)];
            v = (v / p2) * (p2 - 1ULL);
        }
    }

    u64 answer = 0ULL;
    for (u64 m = 1ULL; m <= root; ++m) {
        const u64 m2 = m * m;
        const u64 count = N / m2;
        const u64 add = static_cast<u64>(
            (static_cast<u128>(j2[static_cast<std::size_t>(m)] % kMod) * (count % kMod)) % kMod);
        answer += add;
        if (answer >= kMod) {
            answer -= kMod;
        }
    }

    return answer;
}

}  // namespace

int main() {
    assert(S(10ULL) == 24ULL);
    assert(S(100ULL) == 767ULL);

    std::cout << S(100'000'000'000'000ULL) << '\n';
    return 0;
}
