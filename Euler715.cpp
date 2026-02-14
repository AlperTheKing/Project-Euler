#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <future>
#include <iostream>
#include <limits>
#include <string>
#include <atomic>
#include <thread>
#include <vector>

namespace {

using i64 = long long;
using i128 = __int128_t;

constexpr i64 kMod = 1'000'000'007LL;
constexpr i64 kDefaultN = 1'000'000'000'000LL;
constexpr i64 kCheckpointN1 = 10LL;
constexpr i64 kCheckpointExpected1 = 3'053LL;
constexpr i64 kCheckpointN2 = 100'000LL;
constexpr i64 kCheckpointExpected2 = 157'612'967LL;

struct Options {
    i64 n = kDefaultN;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, i64& value) {
    if (arg.rfind(prefix, 0U) != 0U) return false;
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) return false;

    unsigned long long parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') return false;
        const unsigned long long digit = static_cast<unsigned long long>(ch - '0');
        if (parsed > (std::numeric_limits<unsigned long long>::max() - digit) / 10ULL) return false;
        parsed = parsed * 10ULL + digit;
    }
    if (parsed > static_cast<unsigned long long>(std::numeric_limits<i64>::max())) return false;
    value = static_cast<i64>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }
    return options.n >= 1;
}

inline i64 mod_norm(i64 x) {
    x %= kMod;
    if (x < 0) x += kMod;
    return x;
}

inline i64 mod_add(i64 a, i64 b) {
    i64 s = a + b;
    if (s >= kMod) s -= kMod;
    return s;
}

inline i64 mod_mul(i64 a, i64 b) {
    return static_cast<i64>((static_cast<i128>(a) * b) % kMod);
}

std::vector<int> sieve_primes(const int n) {
    std::vector<char> is_prime(static_cast<std::size_t>(n + 1), 1);
    if (n >= 0) is_prime[0] = 0;
    if (n >= 1) is_prime[1] = 0;
    for (int i = 2; 1LL * i * i <= n; ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) continue;
        for (int j = i * i; j <= n; j += i) is_prime[static_cast<std::size_t>(j)] = 0;
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) primes.push_back(i);
    }
    return primes;
}

std::pair<std::vector<i64>, std::vector<i64>> sum_prime_cubes(
    const i64 n,
    const int L,
    const int XL,
    const std::vector<int>& primes
) {
    std::vector<i64> V(static_cast<std::size_t>(XL), 0);
    i64 s = -1;
    for (int i = 1; i < XL; ++i) {
        const i64 ii = i;
        s = mod_norm(s + mod_mul(mod_mul(ii % kMod, ii % kMod), ii % kMod));
        V[static_cast<std::size_t>(i)] = s;
    }

    std::vector<i64> bigV(static_cast<std::size_t>(L + 1), 0);
    const i64 mod2 = kMod * 2LL;
    for (int i = 1; i <= L; ++i) {
        const i64 x = (n / i) % mod2;
        const i64 y = static_cast<i64>((static_cast<i128>(x) * (x + 1LL) % mod2) / 2);
        bigV[static_cast<std::size_t>(i)] = mod_norm(mod_mul(y % kMod, y % kMod) - 1);
    }

    for (const int p : primes) {
        const i64 sp = V[static_cast<std::size_t>(p - 1)];
        const i64 p_sq = 1LL * p * p;
        const i64 p3 = mod_mul(p % kMod, mod_mul(p % kMod, p % kMod));
        const i64 y = n / p;
        const int iL = static_cast<int>(std::min<i64>(y / p, L));

        for (int i = 1; i <= iL; ++i) {
            const i64 z = y / i;
            const i64 v = (z < XL) ? V[static_cast<std::size_t>(z)] : bigV[static_cast<std::size_t>(i * p)];
            bigV[static_cast<std::size_t>(i)] = mod_norm(
                bigV[static_cast<std::size_t>(i)] - mod_mul(p3, mod_norm(v - sp)));
        }

        for (int x = XL - 1; x > 0; --x) {
            if (x < p_sq) break;
            const int z = x / p;
            V[static_cast<std::size_t>(x)] = mod_norm(
                V[static_cast<std::size_t>(x)] - mod_mul(p3, mod_norm(V[static_cast<std::size_t>(z)] - sp)));
        }
    }

    return {std::move(V), std::move(bigV)};
}

std::pair<std::vector<i64>, std::vector<i64>> prime_balance_mod4(
    const i64 n,
    const int L,
    const int XL,
    const std::vector<int>& primes
) {
    std::vector<i64> V(static_cast<std::size_t>(XL), 0);
    for (int i = 3; i < XL; i += 4) V[static_cast<std::size_t>(i)] = -1;
    for (int i = 4; i < XL; i += 4) V[static_cast<std::size_t>(i)] = -1;

    std::vector<i64> bigV(static_cast<std::size_t>(L + 1), 0);
    for (int i = 1; i <= L; ++i) {
        if (((n / i - 1) % 4) >= 2) bigV[static_cast<std::size_t>(i)] = -1;
    }

    for (std::size_t idx = 1; idx < primes.size(); ++idx) {
        const int p = primes[idx];
        const i64 sp = V[static_cast<std::size_t>(p - 1)];
        const i64 y = n / p;
        const int iL = static_cast<int>(std::min<i64>(y / p, L));
        const i64 f = 2 - (p % 4);

        for (int i = 1; i <= iL; ++i) {
            const i64 z = y / i;
            if (z < XL) {
                bigV[static_cast<std::size_t>(i)] -= f * (V[static_cast<std::size_t>(z)] - sp);
            } else {
                bigV[static_cast<std::size_t>(i)] -= f * (bigV[static_cast<std::size_t>(i * p)] - sp);
            }
        }

        const i64 p_sq = 1LL * p * p;
        for (int x = XL - 1; x > 0; --x) {
            if (x < p_sq) break;
            const int z = x / p;
            V[static_cast<std::size_t>(x)] -= f * (V[static_cast<std::size_t>(z)] - sp);
        }
    }

    return {std::move(V), std::move(bigV)};
}

struct DfsContext {
    const std::vector<int>& primes;
    const std::vector<i64>& p2;
    const std::vector<i64>& p3;
    int L = 0;
    const std::vector<i64>& V;
    const std::vector<i64>& bigV;

    i64 dfs_branch(const int i, const i64 n0, const i64 x0) const {
        const i64 p_sq = p2[static_cast<std::size_t>(i)];
        if (n0 < p_sq) return 0;

        const i64 p = primes[static_cast<std::size_t>(i)];
        const i64 p_cubed = p3[static_cast<std::size_t>(i)];
        i64 n = n0 / p;
        i64 x = x0 * p;
        i64 f = mod_norm(p_cubed + (p % 4) - 2);
        i64 res = 0;

        const i64 pref = (x <= L) ? bigV[static_cast<std::size_t>(x)] : V[static_cast<std::size_t>(n)];
        res = mod_add(res, mod_mul(f, mod_norm(pref - V[static_cast<std::size_t>(p)])));

        while (true) {
            if (n > p_sq) {
                res = mod_add(res, mod_mul(f, dfs(i + 1, n, x)));
            }
            n /= p;
            if (n == 0) break;

            x *= p;
            f = mod_mul(f, p_cubed);
            res = mod_add(res, f);

            if (n > p) {
                const i64 pref2 = (x <= L) ? bigV[static_cast<std::size_t>(x)] : V[static_cast<std::size_t>(n)];
                res = mod_add(res, mod_mul(f, mod_norm(pref2 - V[static_cast<std::size_t>(p)])));
            }
        }
        return res;
    }

    i64 dfs(const int i0, const i64 n0, const i64 x0) const {
        i64 res = 0;
        for (int i = i0; i < static_cast<int>(primes.size()); ++i) {
            if (n0 < p2[static_cast<std::size_t>(i)]) break;
            res = mod_add(res, dfs_branch(i, n0, x0));
        }
        return res;
    }
};

i64 solve(const i64 n) {
    const int L = static_cast<int>(std::sqrt(static_cast<long double>(n) + 0.5L));
    const int XL = static_cast<int>(n / L) + 1;
    const std::vector<int> primes = sieve_primes(L);

    auto fut_prime_balance = std::async(std::launch::async, [&]() {
        return prime_balance_mod4(n, L, XL, primes);
    });
    auto [V1, bigV1] = sum_prime_cubes(n, L, XL, primes);
    auto [V2, bigV2] = fut_prime_balance.get();

    std::vector<i64> V(static_cast<std::size_t>(XL), 0);
    for (int i = 0; i < XL; ++i) {
        V[static_cast<std::size_t>(i)] = mod_norm(
            V1[static_cast<std::size_t>(i)] - mod_norm(V2[static_cast<std::size_t>(i)]));
    }

    std::vector<i64> bigV(static_cast<std::size_t>(L + 1), 0);
    for (int i = 0; i <= L; ++i) {
        bigV[static_cast<std::size_t>(i)] = mod_norm(
            bigV1[static_cast<std::size_t>(i)] - mod_norm(bigV2[static_cast<std::size_t>(i)]));
    }

    std::vector<i64> p2(primes.size(), 0);
    std::vector<i64> p3(primes.size(), 0);
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const i64 p = primes[i];
        p2[i] = p * p;
        p3[i] = mod_mul(p % kMod, mod_mul(p % kMod, p % kMod));
    }

    const DfsContext ctx{primes, p2, p3, L, V, bigV};

    int root_end = 0;
    while (root_end < static_cast<int>(primes.size()) && p2[static_cast<std::size_t>(root_end)] <= n) {
        ++root_end;
    }

    if (root_end == 0) {
        return mod_norm(1 + bigV[1]);
    }

    unsigned int threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    threads = std::min<unsigned int>(threads, static_cast<unsigned int>(root_end));
    if (threads == 1) {
        return mod_norm(1 + bigV[1] + ctx.dfs(0, n, 1));
    }

    std::atomic<int> next_idx(0);
    std::vector<i64> partial(static_cast<std::size_t>(threads), 0);
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(threads));

    for (unsigned int tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid]() {
            i64 local = 0;
            while (true) {
                const int i = next_idx.fetch_add(1, std::memory_order_relaxed);
                if (i >= root_end) break;
                local = mod_add(local, ctx.dfs_branch(i, n, 1));
            }
            partial[static_cast<std::size_t>(tid)] = local;
        });
    }
    for (auto& th : workers) th.join();

    i64 dfs_total = 0;
    for (const i64 v : partial) dfs_total = mod_add(dfs_total, v);
    return mod_norm(1 + bigV[1] + dfs_total);
}

bool run_checkpoints() {
    if (solve(kCheckpointN1) != kCheckpointExpected1) {
        std::cerr << "Validation failed: G(10)\n";
        return false;
    }
    if (solve(kCheckpointN2) != kCheckpointExpected2) {
        std::cerr << "Validation failed: G(10^5)\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }
    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }
    std::cout << solve(options.n) << '\n';
    return 0;
}
