#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int kMaxN = 13;
constexpr int kReachStates = 1 << 16;  // bitset over all suit-subsets (16 subsets)

u64 comb_u64(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }
    u128 num = 1;
    u128 den = 1;
    for (int i = 1; i <= k; ++i) {
        num *= static_cast<u128>(n - k + i);
        den *= static_cast<u128>(i);
    }
    return static_cast<u64>(num / den);
}

std::array<std::array<std::uint16_t, 16>, kReachStates> build_updates() {
    std::array<std::array<std::uint16_t, 16>, kReachStates> update{};

    for (int rmask = 0; rmask < kReachStates; ++rmask) {
        for (int t = 0; t < 16; ++t) {
            int next_mask = rmask;

            for (int matched = 0; matched < 16; ++matched) {
                if ((rmask & (1 << matched)) == 0) {
                    continue;
                }
                const int available = t & (~matched);
                int b = available;
                while (b != 0) {
                    const int suit_bit = __builtin_ctz(b);
                    b &= (b - 1);
                    const int new_subset = matched | (1 << suit_bit);
                    next_mask |= (1 << new_subset);
                }
            }

            update[static_cast<std::size_t>(rmask)][static_cast<std::size_t>(t)] =
                static_cast<std::uint16_t>(next_mask);
        }
    }

    return update;
}

std::array<u64, kMaxN + 1> compute_f_values() {
    const auto update = build_updates();
    std::array<int, 16> popcount{};
    for (int t = 0; t < 16; ++t) {
        popcount[static_cast<std::size_t>(t)] = __builtin_popcount(static_cast<unsigned>(t));
    }

    std::vector<std::vector<u64>> curr(static_cast<std::size_t>(kMaxN + 1),
                                       std::vector<u64>(static_cast<std::size_t>(kReachStates), 0ULL));
    std::vector<std::vector<u64>> next(static_cast<std::size_t>(kMaxN + 1),
                                       std::vector<u64>(static_cast<std::size_t>(kReachStates), 0ULL));

    // Only empty matching set is reachable initially.
    curr[0][1] = 1ULL;  // bit 0 set

    for (int rank = 0; rank < 13; ++rank) {
        for (int n = 0; n <= kMaxN; ++n) {
            std::fill(next[static_cast<std::size_t>(n)].begin(), next[static_cast<std::size_t>(n)].end(), 0ULL);
        }

        for (int n = 0; n <= kMaxN; ++n) {
            for (int rmask = 0; rmask < kReachStates; ++rmask) {
                const u64 ways = curr[static_cast<std::size_t>(n)][static_cast<std::size_t>(rmask)];
                if (ways == 0ULL) {
                    continue;
                }

                for (int t = 0; t < 16; ++t) {
                    const int nn = n + popcount[static_cast<std::size_t>(t)];
                    if (nn > kMaxN) {
                        continue;
                    }
                    const int nr = update[static_cast<std::size_t>(rmask)][static_cast<std::size_t>(t)];
                    next[static_cast<std::size_t>(nn)][static_cast<std::size_t>(nr)] += ways;
                }
            }
        }

        curr.swap(next);
    }

    std::array<u64, kMaxN + 1> f{};
    for (int n = 0; n <= kMaxN; ++n) {
        u64 badugi_hands = 0ULL;
        for (int rmask = 0; rmask < kReachStates; ++rmask) {
            if ((rmask & (1 << 15)) != 0) {  // suit subset {all 4 suits} is reachable
                badugi_hands += curr[static_cast<std::size_t>(n)][static_cast<std::size_t>(rmask)];
            }
        }
        f[static_cast<std::size_t>(n)] = badugi_hands;
    }

    return f;
}

bool run_checkpoints() {
    const auto f = compute_f_values();

    if (f[5] != 514800ULL) {
        std::cerr << "Checkpoint failed: f(5)\n";
        return false;
    }

    // Sanity: for each n <= 13, badugi and non-badugi counts partition all n-card hands.
    // Here we only verify that f(n) never exceeds C(52,n) and f(4)=17160 (all 4-card badugis).
    for (int n = 0; n <= kMaxN; ++n) {
        if (f[static_cast<std::size_t>(n)] > comb_u64(52, n)) {
            std::cerr << "Checkpoint failed: f(n) exceeds total hands for n=" << n << '\n';
            return false;
        }
    }

    if (f[4] != 17160ULL) {
        std::cerr << "Checkpoint failed: f(4)\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    bool skip_checkpoints = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            skip_checkpoints = true;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!skip_checkpoints && !run_checkpoints()) {
        return 2;
    }

    const auto f = compute_f_values();
    u64 answer = 0ULL;
    for (int n = 4; n <= 13; ++n) {
        answer += f[static_cast<std::size_t>(n)];
    }

    std::cout << answer << '\n';
    return 0;
}
