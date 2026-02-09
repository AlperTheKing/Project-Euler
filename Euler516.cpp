#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMask32 = 0xFFFF'FFFFULL;

u64 mod_mul(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * b) % mod);
}

u64 mod_pow(u64 base, u64 exp, const u64 mod) {
    u64 result = 1ULL % mod;
    u64 cur = base % mod;
    u64 e = exp;
    while (e > 0ULL) {
        if (e & 1ULL) {
            result = mod_mul(result, cur, mod);
        }
        cur = mod_mul(cur, cur, mod);
        e >>= 1ULL;
    }
    return result;
}

bool is_prime(const u64 n) {
    if (n < 2ULL) {
        return false;
    }
    for (const u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL}) {
        if (n == p) {
            return true;
        }
        if (n % p == 0ULL) {
            return false;
        }
    }

    u64 d = n - 1ULL;
    int s = 0;
    while ((d & 1ULL) == 0ULL) {
        d >>= 1ULL;
        ++s;
    }

    for (const u64 a : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL}) {
        if (a >= n) {
            continue;
        }
        u64 x = mod_pow(a, d, n);
        if (x == 1ULL || x == n - 1ULL) {
            continue;
        }
        bool witness = true;
        for (int r = 1; r < s; ++r) {
            x = mod_mul(x, x, n);
            if (x == n - 1ULL) {
                witness = false;
                break;
            }
        }
        if (witness) {
            return false;
        }
    }
    return true;
}

std::vector<u64> generate_hamming(const u64 limit) {
    std::vector<u64> out;
    for (u64 a = 1ULL; a <= limit; a *= 2ULL) {
        for (u64 b = a; b <= limit; b *= 3ULL) {
            for (u64 c = b; c <= limit; c *= 5ULL) {
                out.push_back(c);
                if (c > limit / 5ULL) {
                    break;
                }
            }
            if (b > limit / 3ULL) {
                break;
            }
        }
        if (a > limit / 2ULL) {
            break;
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

u64 solve(const u64 limit) {
    const std::vector<u64> hamming = generate_hamming(limit);

    std::vector<u64> prefix_mod(static_cast<std::size_t>(hamming.size() + 1), 0ULL);
    for (std::size_t i = 0; i < hamming.size(); ++i) {
        prefix_mod[i + 1] = (prefix_mod[i] + (hamming[i] & kMask32)) & kMask32;
    }

    std::vector<u64> admissible_primes;
    admissible_primes.reserve(hamming.size());
    for (const u64 h : hamming) {
        const u64 p = h + 1ULL;
        if (p > 5ULL && p <= limit && is_prime(p)) {
            admissible_primes.push_back(p);
        }
    }
    std::sort(admissible_primes.begin(), admissible_primes.end());

    u64 total_mod = 0ULL;

    const auto sum_hamming_mod = [&](const u64 x) -> u64 {
        const auto it = std::upper_bound(hamming.begin(), hamming.end(), x);
        const std::size_t count = static_cast<std::size_t>(it - hamming.begin());
        return prefix_mod[count];
    };

    std::function<void(std::size_t, u64)> dfs = [&](const std::size_t idx, const u64 prod) {
        const u64 h_sum = sum_hamming_mod(limit / prod);
        total_mod = (total_mod + ((prod & kMask32) * h_sum & kMask32)) & kMask32;

        for (std::size_t i = idx; i < admissible_primes.size(); ++i) {
            const u64 p = admissible_primes[i];
            if (prod > limit / p) {
                break;
            }
            dfs(i + 1, prod * p);
        }
    };

    dfs(0, 1ULL);
    return total_mod;
}

u64 brute(const int limit) {
    std::vector<int> phi(static_cast<std::size_t>(limit + 1), 0);
    for (int i = 0; i <= limit; ++i) {
        phi[static_cast<std::size_t>(i)] = i;
    }
    for (int p = 2; p <= limit; ++p) {
        if (phi[static_cast<std::size_t>(p)] != p) {
            continue;
        }
        for (int m = p; m <= limit; m += p) {
            phi[static_cast<std::size_t>(m)] -= phi[static_cast<std::size_t>(m)] / p;
        }
    }

    auto is_smooth = [](int x) -> bool {
        while (x % 2 == 0) {
            x /= 2;
        }
        while (x % 3 == 0) {
            x /= 3;
        }
        while (x % 5 == 0) {
            x /= 5;
        }
        return x == 1;
    };

    u64 sum_mod = 0ULL;
    for (int n = 1; n <= limit; ++n) {
        if (is_smooth(phi[static_cast<std::size_t>(n)])) {
            sum_mod = (sum_mod + static_cast<u64>(n)) & kMask32;
        }
    }
    return sum_mod;
}

bool run_checkpoints() {
    if (solve(100ULL) != 3'728ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    if (solve(5'000ULL) != brute(5'000)) {
        std::cerr << "Checkpoint failed: formula/bruteforce mismatch at 5000\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 limit = 1'000'000'000'000ULL;
    std::cout << solve(limit) << '\n';
    return 0;
}
