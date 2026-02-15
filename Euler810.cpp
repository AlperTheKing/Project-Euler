#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

int bit_length(u64 x) {
    return x == 0 ? 0 : 64 - __builtin_clzll(x);
}

u32 nth_xor_prime_sieve(int M, u64 N) {
    const u64 size = 1ULL << M;
    std::vector<std::uint8_t> is_prime(size, 1);

    u64 count = 1;  // polynomial x
    if (N == 1) return 2u;

    for (u64 i = 3; i < size; i += 2) {
        if (!is_prime[static_cast<std::size_t>(i)]) continue;
        ++count;
        if (count == N) return static_cast<u32>(i);

        const int l = bit_length(i);
        const int j_start = l - 1;
        const int j_end = M - l;
        if (j_start > j_end) continue;

        for (int j = j_start; j <= j_end; ++j) {
            u64 k = i << j;
            const u64 limit = 1ULL << j;
            for (u64 n = 1; n <= limit; ++n) {
                is_prime[static_cast<std::size_t>(k)] = 0;
                k ^= i * (n & (~n + 1ULL));
            }
        }
    }

    return 0u;
}

}  // namespace

int main() {
    const u32 val10 = nth_xor_prime_sieve(10, 10);
    if (val10 != 41u) {
        std::cerr << "Validation failed: " << val10 << "\n";
        return 1;
    }

    std::cout << nth_xor_prime_sieve(27, 5'000'000ULL) << '\n';
    return 0;
}
