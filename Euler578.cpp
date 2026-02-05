#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1) <= x / (r + 1)) ++r;
    while (r > x / r) --r;
    return r;
}

u64 icbrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while (static_cast<u128>(r + 1) * (r + 1) * (r + 1) <= x) ++r;
    while (static_cast<u128>(r) * r * r > x) --r;
    return r;
}

u64 iroot4_u64(u64 x) {
    long double xd = static_cast<long double>(x);
    u64 r = static_cast<u64>(std::sqrt(std::sqrt(xd)));
    while (static_cast<u128>(r + 1) * (r + 1) * (r + 1) * (r + 1) <= x) ++r;
    while (static_cast<u128>(r) * r * r * r > x) --r;
    return r;
}

struct PrimeTable {
    int limit = 0;
    std::vector<int> primes;
    std::vector<int> pi_small;

    explicit PrimeTable(int sieve_limit) { build(sieve_limit); }

    void build(int sieve_limit) {
        limit = sieve_limit;
        std::vector<int> lp(limit + 1, 0);
        primes.reserve(limit / 10);

        for (int i = 2; i <= limit; ++i) {
            if (lp[i] == 0) {
                lp[i] = i;
                primes.push_back(i);
            }
            for (int p : primes) {
                const long long v = 1LL * p * i;
                if (p > lp[i] || v > limit) break;
                lp[static_cast<int>(v)] = p;
            }
        }

        pi_small.assign(limit + 1, 0);
        int ptr = 0;
        for (int i = 1; i <= limit; ++i) {
            pi_small[i] = pi_small[i - 1];
            if (ptr < static_cast<int>(primes.size()) && primes[ptr] == i) {
                ++pi_small[i];
                ++ptr;
            }
        }
    }
};

class PrimeCounter {
public:
    explicit PrimeCounter(const PrimeTable& tbl) : table_(tbl) { build_phi_table(); }

    u64 pi(u64 x) const {
        if (x <= static_cast<u64>(table_.limit)) return static_cast<u64>(table_.pi_small[x]);

        auto& cache = pi_cache();
        auto it = cache.find(x);
        if (it != cache.end()) return it->second;

        const u64 a = pi(iroot4_u64(x));
        const u64 b = pi(isqrt_u64(x));
        const u64 c = pi(icbrt_u64(x));

        i64 sum = static_cast<i64>(phi(x, static_cast<int>(a)));
        sum += static_cast<i64>((b + a - 2) * (b - a + 1) / 2);

        for (u64 i = a + 1; i <= b; ++i) {
            const u64 w = x / static_cast<u64>(table_.primes[static_cast<std::size_t>(i - 1)]);
            sum -= static_cast<i64>(pi(w));
            if (i <= c) {
                const u64 bi = pi(isqrt_u64(w));
                for (u64 j = i; j <= bi; ++j) {
                    const u64 q =
                        w / static_cast<u64>(table_.primes[static_cast<std::size_t>(j - 1)]);
                    sum -= static_cast<i64>(pi(q)) - static_cast<i64>(j - 1);
                }
            }
        }

        const u64 ans = static_cast<u64>(sum);
        cache.emplace(x, ans);
        return ans;
    }

private:
    using i64 = std::int64_t;
    static constexpr int kSmallPhiPrimes = 7;      // 2,3,5,7,11,13,17
    static constexpr int kSmallPhiMod = 510510;    // Product of the above primes.
    static constexpr int kPhiCacheSMax = 128;
    static constexpr u64 kPhiCacheXMax = 50'000'000'000ULL;

    using u32 = std::uint32_t;
    const PrimeTable& table_;
    std::vector<std::array<u32, kSmallPhiPrimes + 1>> phi_table_;

    void build_phi_table() {
        phi_table_.assign(kSmallPhiMod + 1, {});
        for (int n = 0; n <= kSmallPhiMod; ++n) {
            phi_table_[n][0] = static_cast<u32>(n);
        }

        for (int s = 1; s <= kSmallPhiPrimes; ++s) {
            const int p = table_.primes[static_cast<std::size_t>(s - 1)];
            for (int n = 0; n <= kSmallPhiMod; ++n) {
                phi_table_[n][s] = phi_table_[n][s - 1] - phi_table_[n / p][s - 1];
            }
        }
    }

    static std::unordered_map<u64, u64>& pi_cache() {
        static thread_local std::unordered_map<u64, u64> cache;
        return cache;
    }

    static std::unordered_map<u64, u64>& phi_cache() {
        static thread_local std::unordered_map<u64, u64> cache;
        return cache;
    }

    u64 phi(u64 x, int s) const {
        if (s == 0) return x;
        if (s <= kSmallPhiPrimes) {
            const u64 block = x / static_cast<u64>(kSmallPhiMod);
            const u64 rem = x % static_cast<u64>(kSmallPhiMod);
            return block * static_cast<u64>(phi_table_[kSmallPhiMod][s]) +
                   static_cast<u64>(phi_table_[static_cast<std::size_t>(rem)][s]);
        }

        if (x <= static_cast<u64>(table_.limit) &&
            static_cast<u64>(table_.primes[static_cast<std::size_t>(s - 1)]) *
                    static_cast<u64>(table_.primes[static_cast<std::size_t>(s - 1)]) >
                x) {
            return static_cast<u64>(table_.pi_small[static_cast<std::size_t>(x)] - s + 1);
        }

        const u64 ps = static_cast<u64>(table_.primes[static_cast<std::size_t>(s - 1)]);
        if (ps * ps > x) {
            return pi(x) - static_cast<u64>(s) + 1ULL;
        }

        if (s <= kPhiCacheSMax && x <= kPhiCacheXMax) {
            const u64 key = (x << 8) ^ static_cast<u64>(s);
            auto& cache = phi_cache();
            auto it = cache.find(key);
            if (it != cache.end()) return it->second;
            const u64 ans = phi(x, s - 1) - phi(x / ps, s - 1);
            cache.emplace(key, ans);
            return ans;
        }

        return phi(x, s - 1) - phi(x / ps, s - 1);
    }
};

class Solver578 {
public:
    Solver578(const PrimeTable& table, const PrimeCounter& pc) : table_(table), pc_(pc) {}

    u64 solve(u64 n, unsigned threads) const {
        Context base_ctx;
        u64 total = rough_squarefree_count(n, 0, base_ctx);

        std::vector<Task> tasks;
        tasks.reserve(250000);
        build_top_level_tasks(n, tasks);

        if (tasks.empty()) return total;

        if (threads == 0) threads = 1;
        threads = std::min<unsigned>(threads, static_cast<unsigned>(tasks.size()));
        if (threads == 1 || tasks.size() < 2048) {
            Context ctx;
            for (const auto& task : tasks) {
                total += count_from(task.limit, task.start_index, task.max_exp, ctx);
            }
            return total;
        }

        std::atomic<std::size_t> next_task{0};
        std::vector<u64> partial(threads, 0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid]() {
                Context ctx;
                u64 local = 0;
                while (true) {
                    const std::size_t idx = next_task.fetch_add(1, std::memory_order_relaxed);
                    if (idx >= tasks.size()) break;
                    const Task& task = tasks[idx];
                    local += count_from(task.limit, task.start_index, task.max_exp, ctx);
                }
                partial[tid] = local;
            });
        }

        for (auto& th : workers) th.join();
        for (u64 v : partial) total += v;
        return total;
    }

private:
    struct Task {
        u64 limit;
        int start_index;
        int max_exp;
    };

    struct Context {
        std::unordered_map<u64, u64> rough_cache;
    };

    const PrimeTable& table_;
    const PrimeCounter& pc_;

    void build_top_level_tasks(u64 n, std::vector<Task>& tasks) const {
        const auto& primes = table_.primes;
        for (int i = 0; i < static_cast<int>(primes.size()); ++i) {
            const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(i)]);
            if (p > n / p) break;
            u64 pe = p * p;
            for (int e = 2; pe <= n; ++e) {
                tasks.push_back({n / pe, i + 1, e});
                if (pe > n / p) break;
                pe *= p;
            }
        }
    }

    u64 rough_squarefree_count(u64 x, int start_index, Context& ctx) const {
        const u64 m = pc_.pi(x);
        const u64 start_u = static_cast<u64>(start_index);
        if (start_u >= m) return 1;

        const auto& primes = table_.primes;
        const u64 ps = static_cast<u64>(primes[static_cast<std::size_t>(start_index)]);
        if (ps > x / ps) {
            return 1ULL + (m - start_u);
        }

        const u64 key = (x << 20) ^ static_cast<u64>(start_index);
        auto it = ctx.rough_cache.find(key);
        if (it != ctx.rough_cache.end()) return it->second;

        if (ps <= x / ps && ps * ps > x / ps) {
            // At most two primes can appear in the remaining squarefree tail.
            u64 ans = 1ULL + (m - start_u);  // empty set + single primes.
            const u64 r = pc_.pi(isqrt_u64(x));
            for (int i = start_index; i < static_cast<int>(r); ++i) {
                const u64 lim = x / static_cast<u64>(primes[static_cast<std::size_t>(i)]);
                const u64 upto = pc_.pi(lim);
                if (upto > static_cast<u64>(i + 1)) ans += upto - static_cast<u64>(i + 1);
            }
            ctx.rough_cache.emplace(key, ans);
            return ans;
        }

        const u64 r = pc_.pi(isqrt_u64(x));
        if (start_index >= static_cast<int>(r)) {
            const u64 ans = 1ULL + (m - start_u);
            ctx.rough_cache.emplace(key, ans);
            return ans;
        }

        // Recurse only while choosing the next prime keeps room for at least one more factor.
        u64 ans = 1ULL + (m - r);
        for (int i = start_index; i < static_cast<int>(r); ++i) {
            ans += rough_squarefree_count(
                x / static_cast<u64>(primes[static_cast<std::size_t>(i)]), i + 1, ctx);
        }

        ctx.rough_cache.emplace(key, ans);
        return ans;
    }

    u64 count_from(u64 limit, int start_index, int max_exp, Context& ctx) const {
        u64 total = rough_squarefree_count(limit, start_index, ctx);
        const auto& primes = table_.primes;

        for (int i = start_index; i < static_cast<int>(primes.size()); ++i) {
            const u64 p = static_cast<u64>(primes[static_cast<std::size_t>(i)]);
            if (p > limit / p) break;  // p^2 > limit, so exponent >=2 impossible.

            u64 pe = p * p;
            for (int e = 2; e <= max_exp && pe <= limit; ++e) {
                total += count_from(limit / pe, i + 1, e, ctx);
                if (pe > limit / p) break;
                pe *= p;
            }
        }
        return total;
    }
};

bool is_decreasing_prime_power(u64 n) {
    if (n <= 1) return true;
    std::vector<int> exponents;
    u64 x = n;
    for (u64 p = 2; p * p <= x; ++p) {
        if (x % p != 0) continue;
        int cnt = 0;
        while (x % p == 0) {
            x /= p;
            ++cnt;
        }
        exponents.push_back(cnt);
    }
    if (x > 1) exponents.push_back(1);
    for (std::size_t i = 1; i < exponents.size(); ++i) {
        if (exponents[i - 1] < exponents[i]) return false;
    }
    return true;
}

u64 brute_count(u64 n) {
    u64 cnt = 0;
    for (u64 i = 1; i <= n; ++i) {
        if (is_decreasing_prime_power(i)) ++cnt;
    }
    return cnt;
}

bool run_validations(const Solver578& solver, unsigned threads) {
    struct Check {
        u64 n;
        u64 expected;
    };

    const std::vector<Check> checks = {
        {100ULL, 94ULL},
        {1'000'000ULL, 922'052ULL},
    };

    for (const auto& chk : checks) {
        const u64 got = solver.solve(chk.n, threads);
        if (got != chk.expected) {
            std::cerr << "Validation failed for C(" << chk.n << "): got " << got
                      << ", expected " << chk.expected << "\n";
            return false;
        }
    }

    // Extra checkpoint against brute force on a smaller bound.
    const u64 small_n = 10'000ULL;
    const u64 small_expected = brute_count(small_n);
    const u64 small_got = solver.solve(small_n, threads);
    if (small_got != small_expected) {
        std::cerr << "Validation failed for brute-force C(" << small_n << "): got " << small_got
                  << ", expected " << small_expected << "\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    u64 n = 10'000'000'000'000ULL;  // 10^13
    if (argc > 1) n = std::strtoull(argv[1], nullptr, 10);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    if (argc > 2) threads = static_cast<unsigned>(std::strtoul(argv[2], nullptr, 10));
    if (threads == 0) threads = 1;

    const u64 max_needed = std::max<u64>(n, 1'000'000ULL);
    const int sieve_limit =
        static_cast<int>(std::max<u64>(5'000'000ULL, isqrt_u64(max_needed) + 1024ULL));

    PrimeTable table(sieve_limit);
    PrimeCounter prime_counter(table);
    Solver578 solver(table, prime_counter);

    if (!run_validations(solver, threads)) return 1;

    std::cout << solver.solve(n, threads) << "\n";
    return 0;
}
