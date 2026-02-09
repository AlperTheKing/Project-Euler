#include <cstdint>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr i64 kMod = 100000007LL;

i64 mod_pow(i64 base, i64 exp) {
    i64 result = 1 % kMod;
    base %= kMod;
    while (exp > 0) {
        if ((exp & 1LL) != 0LL) {
            result = static_cast<i64>((static_cast<__int128>(result) * base) % kMod);
        }
        base = static_cast<i64>((static_cast<__int128>(base) * base) % kMod);
        exp >>= 1LL;
    }
    return result;
}

i64 comfortable_arrangements_mod(const int n) {
    std::vector<i64> fact(static_cast<std::size_t>(n + 2), 1);
    for (int i = 1; i <= n + 1; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<i64>((static_cast<__int128>(fact[static_cast<std::size_t>(i - 1)]) * i) % kMod);
    }

    std::vector<i64> inv_fact(static_cast<std::size_t>(n + 2), 1);
    inv_fact[static_cast<std::size_t>(n + 1)] = mod_pow(fact[static_cast<std::size_t>(n + 1)], kMod - 2);
    for (int i = n + 1; i >= 1; --i) {
        inv_fact[static_cast<std::size_t>(i - 1)] =
            static_cast<i64>((static_cast<__int128>(inv_fact[static_cast<std::size_t>(i)]) * i) % kMod);
    }

    std::vector<i64> pow2(static_cast<std::size_t>(n + 1), 1);
    for (int i = 1; i <= n; ++i) {
        pow2[static_cast<std::size_t>(i)] = (pow2[static_cast<std::size_t>(i - 1)] * 2LL) % kMod;
    }

    const int boundary_sum[3] = {0, 1, 2};
    const int boundary_mult[3] = {1, 2, 1};

    i64 total = 0;
    for (int m = 1; m <= n; ++m) {
        const i64 fm = fact[static_cast<std::size_t>(m)];
        const i64 fm1 = fact[static_cast<std::size_t>(m - 1)];

        for (int t = 0; t < 3; ++t) {
            const int s = boundary_sum[t];
            const int mult = boundary_mult[t];

            const int g2 = n - (2 * m - 1 + s);
            if (g2 < 0 || g2 > (m - 1)) {
                continue;
            }

            i64 comb = fm1;
            comb = static_cast<i64>((static_cast<__int128>(comb) * inv_fact[static_cast<std::size_t>(g2)]) % kMod);
            comb = static_cast<i64>(
                (static_cast<__int128>(comb) * inv_fact[static_cast<std::size_t>(m - 1 - g2)]) % kMod);

            i64 term = (static_cast<i64>(mult) * comb) % kMod;
            term = static_cast<i64>((static_cast<__int128>(term) * fm) % kMod);
            term = static_cast<i64>((static_cast<__int128>(term) * fm1) % kMod);
            term = static_cast<i64>((static_cast<__int128>(term) * fact[static_cast<std::size_t>(s + g2)]) % kMod);
            term = static_cast<i64>((static_cast<__int128>(term) * pow2[static_cast<std::size_t>(g2)]) % kMod);

            total += term;
            if (total >= kMod) {
                total -= kMod;
            }
        }
    }

    return total;
}

u64 brute_count(const int n) {
    const u64 full_mask = (n == 64) ? ~0ULL : ((1ULL << n) - 1ULL);
    std::unordered_map<u64, u64> memo;
    memo.reserve(1U << n);

    std::function<u64(u64)> dfs = [&](const u64 mask) -> u64 {
        if (mask == full_mask) {
            return 1ULL;
        }
        const auto it = memo.find(mask);
        if (it != memo.end()) {
            return it->second;
        }

        int best = 3;
        std::vector<int> choices;
        for (int i = 0; i < n; ++i) {
            if ((mask & (1ULL << i)) != 0ULL) {
                continue;
            }
            int occ = 0;
            if (i > 0 && (mask & (1ULL << (i - 1))) != 0ULL) {
                ++occ;
            }
            if (i + 1 < n && (mask & (1ULL << (i + 1))) != 0ULL) {
                ++occ;
            }
            if (occ < best) {
                best = occ;
                choices.clear();
                choices.push_back(i);
            } else if (occ == best) {
                choices.push_back(i);
            }
        }

        u64 total = 0;
        for (const int seat : choices) {
            total += dfs(mask | (1ULL << seat));
        }
        memo.emplace(mask, total);
        return total;
    };

    return dfs(0ULL);
}

bool run_checkpoints() {
    if (brute_count(4) != 8ULL) {
        std::cerr << "Checkpoint failed: T(4)\n";
        return false;
    }

    if (comfortable_arrangements_mod(10) != 61632LL) {
        std::cerr << "Checkpoint failed: T(10)\n";
        return false;
    }

    if (comfortable_arrangements_mod(1000) != 47255094LL) {
        std::cerr << "Checkpoint failed: T(1000) mod 100000007\n";
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

    std::cout << comfortable_arrangements_mod(1000000) << '\n';
    return 0;
}
