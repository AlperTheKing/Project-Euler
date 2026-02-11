#include <cassert>
#include <cstdint>
#include <iostream>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 F_bruteforce(const u64 n) {
    u64 best = 0;
    for (u64 m = 0; m <= n; ++m) {
        const u64 v = static_cast<u64>(__builtin_popcountll(m)) +
                      static_cast<u64>(__builtin_popcountll(n - m)) -
                      static_cast<u64>(__builtin_popcountll(n));
        if (v > best) {
            best = v;
        }
    }
    return best;
}

u128 S_fast(const u64 N) {
    if (N == 0U) {
        return 0;
    }

    const int K = 63 - __builtin_clzll(N);
    u128 ans = 0;

    for (int k = 1; k < K; ++k) {
        const u128 p2 = static_cast<u128>(1ULL) << k;
        ans += static_cast<u128>(k - 1) * p2 + 1U;
    }

    const u64 block_start = 1ULL << K;
    const u64 M = N - block_start + 1ULL;
    const u64 A = block_start + 1ULL;
    const u64 B = N + 1ULL;

    u128 sum_min = 0;
    for (int t = 1; t <= K; ++t) {
        const u64 p2 = 1ULL << t;
        const u64 count = B / p2 - (A - 1ULL) / p2;
        sum_min += count;
    }

    ans += static_cast<u128>(K) * static_cast<u128>(M) - sum_min;
    return ans;
}

u64 S_bruteforce(const u64 N) {
    u64 s = 0;
    for (u64 n = 1; n <= N; ++n) {
        s += F_bruteforce(n);
    }
    return s;
}

void print_u128(const u128 x) {
    if (x >= 10U) {
        print_u128(x / 10U);
    }
    const char digit = static_cast<char>('0' + static_cast<int>(x % 10U));
    std::cout << digit;
}

}  // namespace

int main() {
    assert(S_fast(100ULL) == 389ULL);
    assert(S_fast(10'000'000ULL) == 203'222'840ULL);
    assert(S_fast(10'000ULL) == S_bruteforce(10'000ULL));

    const u128 ans = S_fast(10'000'000'000'000'000ULL);
    print_u128(ans);
    std::cout << '\n';
    return 0;
}

