#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;
using u32 = std::uint32_t;

constexpr u64 kMod = 1'000'000'000ULL;

struct PrimeTable {
    int limit = 0;
    std::vector<int> primes;
    std::vector<int> pi;      // pi[x] = #primes <= x
    std::vector<u64> ps_mod;  // sum_{p<=x} p mod kMod
};

PrimeTable sieve_primes(const int limit) {
    PrimeTable out;
    out.limit = limit;
    std::vector<bool> is_prime(static_cast<std::size_t>(limit + 1), true);
    if (limit >= 0) {
        is_prime[0] = false;
    }
    if (limit >= 1) {
        is_prime[1] = false;
    }
    for (int p = 2; (i128)p * p <= limit; ++p) {
        if (!is_prime[static_cast<std::size_t>(p)]) {
            continue;
        }
        for (int m = p * p; m <= limit; m += p) {
            is_prime[static_cast<std::size_t>(m)] = false;
        }
    }
    out.pi.assign(static_cast<std::size_t>(limit + 1), 0);
    out.ps_mod.assign(static_cast<std::size_t>(limit + 1), 0ULL);
    int cnt = 0;
    u64 sum_mod = 0ULL;
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            out.primes.push_back(i);
            ++cnt;
            sum_mod += static_cast<u64>(i);
            sum_mod %= kMod;
        }
        out.pi[static_cast<std::size_t>(i)] = cnt;
        out.ps_mod[static_cast<std::size_t>(i)] = sum_mod;
    }
    return out;
}

struct Min25PrimeSummatory {
    u64 n = 0;
    u64 sqrt_n = 0;
    std::vector<u64> vals;   // distinct n / i
    std::vector<u64> g_cnt;  // pi(v)
    std::vector<u64> g_sum;  // sum_{p<=v} p mod kMod
    std::vector<int> id_small;
    std::vector<int> id_large;

    int id_of(const u64 x) const {
        if (x <= sqrt_n) {
            return id_small[static_cast<std::size_t>(x)];
        }
        return id_large[static_cast<std::size_t>(n / x)];
    }

    u64 pi_of(const u64 x) const { return g_cnt[static_cast<std::size_t>(id_of(x))]; }
    u64 prime_sum_mod(const u64 x) const { return g_sum[static_cast<std::size_t>(id_of(x))]; }
};

Min25PrimeSummatory min25_prime_pi_sum(const u64 n, const std::vector<int>& primes) {
    Min25PrimeSummatory ms;
    ms.n = n;
    ms.sqrt_n = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((ms.sqrt_n + 1) * (ms.sqrt_n + 1) <= n) {
        ++ms.sqrt_n;
    }
    while (ms.sqrt_n * ms.sqrt_n > n) {
        --ms.sqrt_n;
    }

    for (u64 l = 1ULL; l <= n;) {
        const u64 v = n / l;
        const u64 r = n / v;
        ms.vals.push_back(v);
        l = r + 1ULL;
    }

    const int m = static_cast<int>(ms.vals.size());
    ms.g_cnt.assign(static_cast<std::size_t>(m), 0ULL);
    ms.g_sum.assign(static_cast<std::size_t>(m), 0ULL);
    ms.id_small.assign(static_cast<std::size_t>(ms.sqrt_n + 1), -1);
    ms.id_large.assign(static_cast<std::size_t>(ms.sqrt_n + 1), -1);

    for (int i = 0; i < m; ++i) {
        const u64 v = ms.vals[static_cast<std::size_t>(i)];
        if (v <= ms.sqrt_n) {
            ms.id_small[static_cast<std::size_t>(v)] = i;
        } else {
            ms.id_large[static_cast<std::size_t>(n / v)] = i;
        }

        ms.g_cnt[static_cast<std::size_t>(i)] = v - 1ULL;  // count of integers in [2..v]
        const u128 vv = static_cast<u128>(v);
        const u128 sum_1_to_v = vv * (vv + 1U) / 2U;
        const u64 sum_mod = static_cast<u64>(sum_1_to_v % static_cast<u128>(kMod));
        ms.g_sum[static_cast<std::size_t>(i)] = (sum_mod + kMod - 1ULL) % kMod;  // sum_{2..v} i
    }

    u64 primes_count = 0ULL;   // pi(p-1)
    u64 primes_sum_mod = 0ULL;  // sum_{q<p} q mod kMod

    for (const int p : primes) {
        const u64 pu = static_cast<u64>(p);
        const u64 p2 = pu * pu;
        if (p2 > n) {
            break;
        }

        for (int i = 0; i < m; ++i) {
            const u64 v = ms.vals[static_cast<std::size_t>(i)];
            if (v < p2) {
                break;
            }
            const u64 vp = v / pu;
            const int j = ms.id_of(vp);

            const u64 cnt_j = ms.g_cnt[static_cast<std::size_t>(j)];
            ms.g_cnt[static_cast<std::size_t>(i)] -= (cnt_j - primes_count);

            const u64 sum_j = ms.g_sum[static_cast<std::size_t>(j)];
            const u64 tmp = (sum_j >= primes_sum_mod) ? (sum_j - primes_sum_mod)
                                                      : (sum_j + kMod - primes_sum_mod);
            const u64 delta = (pu % kMod) * tmp % kMod;
            const u64 cur = ms.g_sum[static_cast<std::size_t>(i)];
            ms.g_sum[static_cast<std::size_t>(i)] = (cur >= delta) ? (cur - delta) : (cur + kMod - delta);
        }

        ++primes_count;
        primes_sum_mod += pu % kMod;
        primes_sum_mod %= kMod;
    }

    return ms;
}

struct PhiComputer {
    static constexpr int kWheelA = 6;
    static constexpr int kWheelMod = 2 * 3 * 5 * 7 * 11 * 13;  // 30030

    const PrimeTable& small;
    const Min25PrimeSummatory& ms;

    std::vector<int> wheel_prefix;  // prefix of residues coprime to wheel modulus
    std::unordered_map<u64, u64> memo;

    explicit PhiComputer(const PrimeTable& small_, const Min25PrimeSummatory& ms_)
        : small(small_), ms(ms_) {
        wheel_prefix.assign(kWheelMod + 1, 0);
        int cnt = 0;
        for (int i = 1; i <= kWheelMod; ++i) {
            if (std::gcd(i, kWheelMod) == 1) {
                ++cnt;
            }
            wheel_prefix[i] = cnt;
        }
        memo.reserve(1 << 20);
    }

    u64 pi_of(const u64 x) const {
        if (x <= static_cast<u64>(small.limit)) {
            return static_cast<u64>(small.pi[static_cast<std::size_t>(x)]);
        }
        return ms.pi_of(x);
    }

    u64 phi_wheel_6(const u64 x) const {
        if (x == 0ULL) {
            return 0ULL;
        }
        const u64 full = x / static_cast<u64>(kWheelMod);
        const int rem = static_cast<int>(x % static_cast<u64>(kWheelMod));
        return full * static_cast<u64>(wheel_prefix[kWheelMod]) + static_cast<u64>(wheel_prefix[rem]);
    }

    u64 phi(u64 x, int a) {
        if (x == 0ULL) {
            return 0ULL;
        }
        if (a <= 0) {
            return x;
        }
        if (a == kWheelA) {
            return phi_wheel_6(x);
        }
        if (a < kWheelA) {
            return phi(x, a - 1) - phi(x / static_cast<u64>(small.primes[static_cast<std::size_t>(a - 1)]),
                                       a - 1);
        }

        const u64 pix = pi_of(x);
        if (static_cast<u64>(a) >= pix) {
            return 1ULL;
        }
        const u64 pa = static_cast<u64>(small.primes[static_cast<std::size_t>(a - 1)]);
        if (pa * pa > x) {
            // Only 1 and primes > pa remain.
            return 1ULL + (pix - static_cast<u64>(a));
        }

        const u64 key = (x << 11) | static_cast<u64>(a);
        const auto it = memo.find(key);
        if (it != memo.end()) {
            return it->second;
        }

        const u64 result = phi(x, a - 1) - phi(x / pa, a - 1);
        memo.emplace(key, result);
        return result;
    }
};

u64 solve(const u64 n) {
    const u64 limit = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    const PrimeTable small = sieve_primes(static_cast<int>(limit));
    const Min25PrimeSummatory ms = min25_prime_pi_sum(n, small.primes);

    // cbrt(n)
    u64 c = static_cast<u64>(std::cbrt(static_cast<long double>(n)));
    while ((c + 1) * (c + 1) * (c + 1) <= n) {
        ++c;
    }
    while (c * c * c > n) {
        --c;
    }

    PhiComputer phi(small, ms);

    u64 ans = ms.prime_sum_mod(n);
    const u64 sum_to_limit = small.ps_mod[static_cast<std::size_t>(limit)];
    ans = (ans + kMod - sum_to_limit) % kMod;  // primes in (limit..n]

    for (std::size_t idx = 0; idx < small.primes.size(); ++idx) {
        const u64 p = static_cast<u64>(small.primes[idx]);
        if (p > limit) {
            break;
        }
        const u64 x = n / p;
        u64 cnt = 0ULL;
        if (p > c) {
            // phi(x, idx) shortcut since p^2 > x (equivalently p^3 > n).
            const u64 pi_x = ms.pi_of(x);
            cnt = 1ULL + (pi_x - static_cast<u64>(idx));
        } else {
            cnt = phi.phi(x, static_cast<int>(idx));
        }
        const u64 contrib = (p % kMod) * (cnt % kMod) % kMod;
        ans += contrib;
        ans %= kMod;
    }

    return ans;
}

u64 brute_small(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n + 1), 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = i;
            if ((i128)i * i <= n) {
                for (int m = i * i; m <= n; m += i) {
                    if (spf[static_cast<std::size_t>(m)] == 0) {
                        spf[static_cast<std::size_t>(m)] = i;
                    }
                }
            }
        }
    }
    u64 sum = 0;
    for (int i = 2; i <= n; ++i) {
        sum += static_cast<u64>(spf[static_cast<std::size_t>(i)]);
    }
    return sum;
}

bool run_checkpoints() {
    if (solve(100ULL) != 1'257ULL) {
        std::cerr << "Checkpoint failed: S(100)\n";
        return false;
    }
    if (solve(100ULL) != brute_small(100)) {
        std::cerr << "Checkpoint failed: solve/brute mismatch at 100\n";
        return false;
    }
    return true;
}

}  // namespace

int main() {
    if (!run_checkpoints()) {
        return 1;
    }

    constexpr u64 n = 1'000'000'000'000ULL;
    const u64 answer = solve(n) % kMod;
    std::cout << std::setw(9) << std::setfill('0') << answer << '\n';
    return 0;
}

