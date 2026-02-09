#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'234'567'890ULL;

int v2_exponent(u64 x) {
    int e = 0;
    while ((x & 1ULL) == 0ULL) {
        x >>= 1ULL;
        ++e;
    }
    return e;
}

bool verify_grundy_pattern(const int n_max) {
    std::vector<int> sg(static_cast<std::size_t>(n_max + 1), 0);
    for (int n = 1; n <= n_max; ++n) {
        std::array<bool, 32> seen{};
        for (int d = 1; d < n; ++d) {
            if (n % d != 0) {
                continue;
            }
            const int g = sg[static_cast<std::size_t>(n - d)];
            if (g < static_cast<int>(seen.size())) {
                seen[static_cast<std::size_t>(g)] = true;
            }
        }
        int mex = 0;
        while (mex < static_cast<int>(seen.size()) && seen[static_cast<std::size_t>(mex)]) {
            ++mex;
        }
        sg[static_cast<std::size_t>(n)] = mex;
        if (mex != v2_exponent(static_cast<u64>(n))) {
            std::cerr << "Grundy mismatch at n=" << n << '\n';
            return false;
        }
    }
    return true;
}

std::array<u64, 64> valuation_counts(const u64 n) {
    std::array<u64, 64> count{};
    for (int k = 0; k < 63; ++k) {
        const u64 a = n >> k;
        if (a == 0ULL) {
            break;
        }
        const u64 b = n >> (k + 1);
        count[static_cast<std::size_t>(k)] = a - b;
    }
    return count;
}

u64 winning_positions_mod(const u64 n, const u64 mod) {
    const auto count = valuation_counts(n);

    u64 losing_mod = 0ULL;
    for (int i = 0; i < 64; ++i) {
        if (count[static_cast<std::size_t>(i)] == 0ULL) {
            continue;
        }
        for (int j = 0; j < 64; ++j) {
            if (count[static_cast<std::size_t>(j)] == 0ULL) {
                continue;
            }
            const int k = i ^ j;
            if (k >= 64 || count[static_cast<std::size_t>(k)] == 0ULL) {
                continue;
            }
            const u64 term = static_cast<u64>(
                (static_cast<u128>(count[static_cast<std::size_t>(i)] % mod) *
                 (count[static_cast<std::size_t>(j)] % mod) % mod) *
                (count[static_cast<std::size_t>(k)] % mod) % mod);
            losing_mod += term;
            losing_mod %= mod;
        }
    }

    const u64 n_mod = n % mod;
    const u64 total_mod =
        static_cast<u64>((static_cast<u128>(n_mod) * n_mod % mod) * n_mod % mod);
    return (total_mod + mod - losing_mod) % mod;
}

u64 winning_positions_exact_small(const u64 n) {
    const auto count = valuation_counts(n);
    u128 losing = 0;
    for (int i = 0; i < 64; ++i) {
        for (int j = 0; j < 64; ++j) {
            const int k = i ^ j;
            if (k >= 64) {
                continue;
            }
            losing += static_cast<u128>(count[static_cast<std::size_t>(i)]) *
                      count[static_cast<std::size_t>(j)] * count[static_cast<std::size_t>(k)];
        }
    }
    const u128 total = static_cast<u128>(n) * n * n;
    return static_cast<u64>(total - losing);
}

u64 brute_winning_positions(const int n) {
    std::vector<int> sg(static_cast<std::size_t>(n + 1), 0);
    for (int x = 1; x <= n; ++x) {
        std::array<bool, 32> seen{};
        for (int d = 1; d < x; ++d) {
            if (x % d == 0) {
                const int g = sg[static_cast<std::size_t>(x - d)];
                if (g < static_cast<int>(seen.size())) {
                    seen[static_cast<std::size_t>(g)] = true;
                }
            }
        }
        int mex = 0;
        while (mex < static_cast<int>(seen.size()) && seen[static_cast<std::size_t>(mex)]) {
            ++mex;
        }
        sg[static_cast<std::size_t>(x)] = mex;
    }

    u64 wins = 0ULL;
    for (int a = 1; a <= n; ++a) {
        for (int b = 1; b <= n; ++b) {
            for (int c = 1; c <= n; ++c) {
                if ((sg[static_cast<std::size_t>(a)] ^ sg[static_cast<std::size_t>(b)] ^
                     sg[static_cast<std::size_t>(c)]) != 0) {
                    ++wins;
                }
            }
        }
    }
    return wins;
}

bool run_checkpoints() {
    if (!verify_grundy_pattern(256)) {
        return false;
    }
    if (winning_positions_exact_small(10ULL) != 692ULL) {
        std::cerr << "Checkpoint failed: S(10)\n";
        return false;
    }
    if (winning_positions_exact_small(100ULL) != 735'494ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    if (winning_positions_exact_small(20ULL) != brute_winning_positions(20)) {
        std::cerr << "Checkpoint failed: exact/bruteforce mismatch at n=20\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 123'456'787'654'321ULL;
    std::cout << winning_positions_mod(n, kMod) << '\n';
    return 0;
}
