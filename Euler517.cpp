#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using u32 = std::uint32_t;
using i64 = std::int64_t;

constexpr u64 kMod = 1'000'000'007ULL;

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = static_cast<u64>((static_cast<u128>(result) * cur) % mod);
        }
        cur = static_cast<u64>((static_cast<u128>(cur) * cur) % mod);
        e >>= 1ULL;
    }
    return result;
}

std::vector<int> small_primes_up_to(const int n) {
    std::vector<bool> is_prime(static_cast<std::size_t>(n + 1), true);
    if (n >= 0) {
        is_prime[0] = false;
    }
    if (n >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; p * p <= n; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int m = p * p; m <= n; m += p) {
            is_prime[static_cast<std::size_t>(m)] = false;
        }
    }
    std::vector<int> out;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            out.push_back(i);
        }
    }
    return out;
}

std::vector<int> segmented_primes(const int low, const int high_exclusive) {
    if (high_exclusive <= low) {
        return {};
    }

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(high_exclusive - 1))) + 1;
    const std::vector<int> base_primes = small_primes_up_to(root);

    const int size = high_exclusive - low;
    std::vector<bool> is_prime(static_cast<std::size_t>(size), true);

    for (const int p : base_primes) {
        i64 start = static_cast<i64>((low + p - 1) / p) * p;
        const i64 pp = static_cast<i64>(p) * p;
        if (start < pp) {
            start = pp;
        }
        for (i64 x = start; x < high_exclusive; x += p) {
            is_prime[static_cast<std::size_t>(x - low)] = false;
        }
    }

    if (low == 0) {
        is_prime[0] = false;
        if (size > 1) {
            is_prime[1] = false;
        }
    } else if (low == 1) {
        is_prime[0] = false;
    }

    std::vector<int> primes;
    for (int i = 0; i < size; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(low + i);
        }
    }
    return primes;
}

u64 ceil_k_sqrt_n(const int n, const int k, const long double root) {
    long double approx = std::floor(static_cast<long double>(k) * root);
    u64 x = static_cast<u64>(approx);
    const u128 target = static_cast<u128>(k) * static_cast<u128>(k) * static_cast<u128>(n);

    while (static_cast<u128>(x + 1ULL) * static_cast<u128>(x + 1ULL) <= target) {
        ++x;
    }
    while (static_cast<u128>(x) * static_cast<u128>(x) > target) {
        --x;
    }
    if (static_cast<u128>(x) * static_cast<u128>(x) == target) {
        return x;
    }
    return x + 1ULL;
}

u64 nCk_mod(const int n, const int k, const std::vector<u32>& fact, const std::vector<u32>& invfact) {
    if (k < 0 || k > n) {
        return 0ULL;
    }
    return static_cast<u64>(fact[static_cast<std::size_t>(n)]) *
           invfact[static_cast<std::size_t>(k)] % kMod *
           invfact[static_cast<std::size_t>(n - k)] % kMod;
}

u64 G_mod(const int n, const std::vector<u32>& fact, const std::vector<u32>& invfact) {
    const long double root = std::sqrt(static_cast<long double>(n));
    const int k_max = static_cast<int>(std::floor(root));

    u64 total = 0ULL;
    for (int k = 0; k <= k_max; ++k) {
        const u64 ceil_term = ceil_k_sqrt_n(n, k, root);
        const int m = static_cast<int>(static_cast<u64>(n) - ceil_term + static_cast<u64>(k));
        if (m < k) {
            continue;
        }
        total += nCk_mod(m, k, fact, invfact);
        total %= kMod;
    }
    return total;
}

u64 solve_sum_primes(const int low_exclusive, const int high_exclusive) {
    const int max_n = high_exclusive - 1;
    std::vector<u32> fact(static_cast<std::size_t>(max_n + 1), 1U);
    std::vector<u32> invfact(static_cast<std::size_t>(max_n + 1), 1U);

    for (int i = 1; i <= max_n; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<u32>(static_cast<u64>(fact[static_cast<std::size_t>(i - 1)]) * i % kMod);
    }
    invfact[static_cast<std::size_t>(max_n)] =
        static_cast<u32>(mod_pow(fact[static_cast<std::size_t>(max_n)], kMod - 2ULL, kMod));
    for (int i = max_n; i >= 1; --i) {
        invfact[static_cast<std::size_t>(i - 1)] =
            static_cast<u32>(static_cast<u64>(invfact[static_cast<std::size_t>(i)]) * i % kMod);
    }

    const std::vector<int> primes = segmented_primes(low_exclusive + 1, high_exclusive);
    u64 total = 0ULL;
    for (const int p : primes) {
        total += G_mod(p, fact, invfact);
        total %= kMod;
    }
    return total;
}

u64 brute_G_small(const int n) {
    const long double a = std::sqrt(static_cast<long double>(n));
    std::unordered_map<u64, u64> memo;
    memo.reserve(2048);

    std::function<u64(int, int)> dfs = [&](const int i, const int j) -> u64 {
        const long double x =
            static_cast<long double>(n) - static_cast<long double>(i) - static_cast<long double>(j) * a;
        if (x < a - 1e-15L) {
            return 1ULL;
        }
        const u64 key = (static_cast<u64>(static_cast<u32>(i)) << 32) |
                        static_cast<u64>(static_cast<u32>(j));
        const auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }
        const u64 value = dfs(i + 1, j) + dfs(i, j + 1);
        memo.emplace(key, value);
        return value;
    };

    return dfs(0, 0);
}

bool run_checkpoints() {
    constexpr int checkpoint_n = 90;
    std::vector<u32> fact(static_cast<std::size_t>(checkpoint_n + 1), 1U);
    std::vector<u32> invfact(static_cast<std::size_t>(checkpoint_n + 1), 1U);
    for (int i = 1; i <= checkpoint_n; ++i) {
        fact[static_cast<std::size_t>(i)] =
            static_cast<u32>(static_cast<u64>(fact[static_cast<std::size_t>(i - 1)]) * i % kMod);
    }
    invfact[static_cast<std::size_t>(checkpoint_n)] =
        static_cast<u32>(mod_pow(fact[static_cast<std::size_t>(checkpoint_n)], kMod - 2ULL, kMod));
    for (int i = checkpoint_n; i >= 1; --i) {
        invfact[static_cast<std::size_t>(i - 1)] =
            static_cast<u32>(static_cast<u64>(invfact[static_cast<std::size_t>(i)]) * i % kMod);
    }

    if (G_mod(90, fact, invfact) != 7'564'511ULL) {
        std::cerr << "Checkpoint failed: G(90)\n";
        return false;
    }
    for (int n = 2; n <= 30; ++n) {
        if (G_mod(n, fact, invfact) != brute_G_small(n) % kMod) {
            std::cerr << "Checkpoint failed: formula/bruteforce mismatch at n=" << n << '\n';
            return false;
        }
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr int low = 10'000'000;
    constexpr int high = 10'010'000;
    std::cout << solve_sum_primes(low, high) << '\n';
    return 0;
}
