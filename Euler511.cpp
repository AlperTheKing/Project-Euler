#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'000'000'000ULL;

std::vector<std::pair<u64, int>> factorize(u64 n) {
    std::vector<std::pair<u64, int>> factors;
    u64 x = n;
    for (u64 p = 2ULL; p * p <= x; p += (p == 2ULL ? 1ULL : 2ULL)) {
        if (x % p != 0ULL) {
            continue;
        }
        int e = 0;
        while (x % p == 0ULL) {
            x /= p;
            ++e;
        }
        factors.push_back({p, e});
    }
    if (x > 1ULL) {
        factors.push_back({x, 1});
    }
    return factors;
}

void build_divisors_rec(const std::vector<std::pair<u64, int>>& factors, const int idx,
                        const u64 cur, std::vector<u64>& out) {
    if (idx == static_cast<int>(factors.size())) {
        out.push_back(cur);
        return;
    }
    const auto [p, e] = factors[static_cast<std::size_t>(idx)];
    u64 value = 1ULL;
    for (int i = 0; i <= e; ++i) {
        build_divisors_rec(factors, idx + 1, cur * value, out);
        value *= p;
    }
}

std::vector<u64> divisors_of(u64 n) {
    const auto factors = factorize(n);
    std::vector<u64> divisors;
    build_divisors_rec(factors, 0, 1ULL, divisors);
    return divisors;
}

std::vector<u32> cyclic_convolution(const std::vector<u32>& a, const std::vector<u32>& b,
                                    const u64 mod) {
    const int k = static_cast<int>(a.size());
    std::vector<u64> accum(static_cast<std::size_t>(k), 0ULL);
    constexpr u64 kReduceThreshold = (1ULL << 62);

    for (int i = 0; i < k; ++i) {
        const u64 ai = a[static_cast<std::size_t>(i)];
        if (ai == 0ULL) {
            continue;
        }
        const int split = k - i;
        for (int j = 0; j < split; ++j) {
            u64& cell = accum[static_cast<std::size_t>(i + j)];
            cell += ai * static_cast<u64>(b[static_cast<std::size_t>(j)]);
            if (cell >= kReduceThreshold) {
                cell %= mod;
            }
        }
        for (int j = split; j < k; ++j) {
            u64& cell = accum[static_cast<std::size_t>(i + j - k)];
            cell += ai * static_cast<u64>(b[static_cast<std::size_t>(j)]);
            if (cell >= kReduceThreshold) {
                cell %= mod;
            }
        }
    }

    std::vector<u32> out(static_cast<std::size_t>(k), 0U);
    for (int i = 0; i < k; ++i) {
        out[static_cast<std::size_t>(i)] =
            static_cast<u32>(accum[static_cast<std::size_t>(i)] % mod);
    }
    return out;
}

std::vector<u32> cyclic_power(std::vector<u32> base, u64 exp, const u64 mod) {
    const int k = static_cast<int>(base.size());
    std::vector<u32> result(static_cast<std::size_t>(k), 0U);
    result[0] = 1U;

    while (exp > 0ULL) {
        if (exp & 1ULL) {
            result = cyclic_convolution(result, base, mod);
        }
        exp >>= 1ULL;
        if (exp > 0ULL) {
            base = cyclic_convolution(base, base, mod);
        }
    }
    return result;
}

u64 seq_mod(const u64 n, const int k) {
    const std::vector<u64> divisors = divisors_of(n);
    std::vector<u32> base(static_cast<std::size_t>(k), 0U);
    for (const u64 d : divisors) {
        const int residue = static_cast<int>(d % static_cast<u64>(k));
        base[static_cast<std::size_t>(residue)] += 1U;
    }

    const std::vector<u32> ways = cyclic_power(base, n, kMod);
    const int target = static_cast<int>((static_cast<u64>(k) - (n % static_cast<u64>(k))) %
                                        static_cast<u64>(k));
    return static_cast<u64>(ways[static_cast<std::size_t>(target)]);
}

u64 seq_mod_small_dp(const int n, const int k) {
    const std::vector<u64> divisors = divisors_of(static_cast<u64>(n));
    std::vector<u64> dp(static_cast<std::size_t>(k), 0ULL);
    dp[0] = 1ULL;

    for (int step = 0; step < n; ++step) {
        std::vector<u64> next(static_cast<std::size_t>(k), 0ULL);
        for (int r = 0; r < k; ++r) {
            if (dp[static_cast<std::size_t>(r)] == 0ULL) {
                continue;
            }
            for (const u64 d : divisors) {
                const int nr = (r + static_cast<int>(d % static_cast<u64>(k))) % k;
                next[static_cast<std::size_t>(nr)] += dp[static_cast<std::size_t>(r)];
            }
        }
        dp.swap(next);
    }
    const int target = (k - (n % k)) % k;
    return dp[static_cast<std::size_t>(target)] % kMod;
}

bool run_checkpoints() {
    if (seq_mod(3ULL, 4) != 4ULL) {
        std::cerr << "Checkpoint failed: Seq(3,4)\n";
        return false;
    }
    if (seq_mod(4ULL, 11) != 8ULL) {
        std::cerr << "Checkpoint failed: Seq(4,11)\n";
        return false;
    }
    if (seq_mod(1'111ULL, 24) != 840'643'584ULL) {
        std::cerr << "Checkpoint failed: Seq(1111,24) mod 1e9\n";
        return false;
    }
    if (seq_mod(8ULL, 15) != seq_mod_small_dp(8, 15)) {
        std::cerr << "Checkpoint failed: fast/DP mismatch\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 1'234'567'898'765ULL;
    constexpr int k = 4'321;
    const u64 answer = seq_mod(n, k);
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}
