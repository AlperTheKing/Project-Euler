#include <cstdint>
#include <iomanip>
#include <iostream>
#include <unordered_map>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

u64 count_configurations(int n, int D) {
    // Place exactly one queen per row. A queen attacks vertically/diagonally only up to D rows away.
    // Row-by-row DP with state = last D columns (4 bits each, most recent at LSB).
    const u64 all_cols = (n == 64) ? ~0ULL : ((1ULL << n) - 1ULL);
    const u64 mask = (D == 0) ? 0ULL : ((1ULL << (4 * D)) - 1ULL);

    std::unordered_map<u64, u64> dp, next;
    dp.reserve(1 << 12);
    next.reserve(1 << 12);
    dp.emplace(0ULL, 1ULL);

    for (int row = 0; row < n; ++row) {
        next.clear();
        const int len = (row < D) ? row : D;
        for (const auto& [state, ways] : dp) {
            u64 blocked = 0;
            for (int dist = 1; dist <= len; ++dist) {
                const int c_prev = static_cast<int>((state >> (4 * (dist - 1))) & 0xFULL);
                blocked |= 1ULL << c_prev;
                const int c1 = c_prev + dist;
                const int c2 = c_prev - dist;
                if (c1 >= 0 && c1 < n) {
                    blocked |= 1ULL << c1;
                }
                if (c2 >= 0 && c2 < n) {
                    blocked |= 1ULL << c2;
                }
            }
            u64 avail = all_cols & ~blocked;
            while (avail) {
                const int col = __builtin_ctzll(avail);
                avail &= (avail - 1);
                const u64 ns = (((state << 4) | static_cast<u64>(col)) & mask);
                next[ns] += ways;
            }
        }
        dp.swap(next);
    }

    u128 total = 0;
    for (const auto& [state, ways] : dp) {
        (void)state;
        total += static_cast<u128>(ways);
    }
    return static_cast<u64>(total);
}

u64 S(int n) {
    u128 sum = 0;
    for (int w = 0; w < n; ++w) {
        const int D = (n - 1) - w;
        sum += static_cast<u128>(count_configurations(n, D));
    }
    return static_cast<u64>(sum);
}

bool run_checkpoints() {
    if (count_configurations(4, 3) != 2ULL) {  // w=0 => D=3
        std::cerr << "Checkpoint failed: Q(4,0)\n";
        return false;
    }
    if (count_configurations(4, 1) != 16ULL) {  // w=2 => D=1
        std::cerr << "Checkpoint failed: Q(4,2)\n";
        return false;
    }
    if (count_configurations(4, 0) != 256ULL) {  // w=3 => D=0
        std::cerr << "Checkpoint failed: Q(4,3)\n";
        return false;
    }
    if (S(4) != 276ULL) {
        std::cerr << "Checkpoint failed: S(4)\n";
        return false;
    }
    if (S(5) != 3347ULL) {
        std::cerr << "Checkpoint failed: S(5)\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }
    std::cout << S(14) << '\n';
    return 0;
}

