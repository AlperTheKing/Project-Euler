#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

namespace {

using u64 = uint64_t;
using u128 = __uint128_t;

constexpr u64 kTargetN = 10000000000000000ULL;
constexpr u64 kStep = 2520ULL;

// Strong lower bound found by a deep top-window scan.
constexpr u64 kSeedBestN = 9999452799100481ULL;
constexpr u64 kSeedBestG = 49580620131241271ULL;

constexpr long double kCTop = 12503.0L / 2520.0L;
constexpr long double kK311 = 50327.0L / 2520.0L;
constexpr long double kK2201 = 49697.0L / 2520.0L;

u64 mod_mul(u64 a, u64 b, u64 mod) { return static_cast<u64>((static_cast<u128>(a) * b) % mod); }

u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 r = 1 % mod;
    while (e > 0) {
        if (e & 1ULL) r = mod_mul(r, a, mod);
        a = mod_mul(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

bool is_prime(u64 n) {
    if (n < 2) return false;

    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n % p == 0) return n == p;
    }

    u64 d = n - 1;
    u64 s = 0;
    while ((d & 1ULL) == 0) {
        d >>= 1ULL;
        ++s;
    }

    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        u64 x = mod_pow(a, d, n);
        if (x == 1 || x == n - 1) continue;

        bool composite = true;
        for (u64 r = 1; r < s; ++r) {
            x = mod_mul(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;
    }

    return true;
}

thread_local std::mt19937_64 g_rng(std::random_device{}());

u64 pollard_rho(u64 n) {
    if ((n & 1ULL) == 0) return 2ULL;
    if (n % 3ULL == 0) return 3ULL;

    while (true) {
        const u64 c = g_rng() % (n - 1) + 1;
        u64 x = g_rng() % (n - 2) + 2;
        u64 y = x;
        u64 d = 1;

        auto f = [&](u64 v) { return (mod_mul(v, v, n) + c) % n; };

        while (d == 1) {
            x = f(x);
            y = f(f(y));
            const u64 diff = (x > y) ? (x - y) : (y - x);
            d = std::gcd(diff, n);
        }

        if (d != n) return d;
    }
}

void factor_rec(u64 n, std::vector<u64>& factors) {
    if (n == 1) return;
    if (is_prime(n)) {
        factors.push_back(n);
        return;
    }

    const u64 d = pollard_rho(n);
    factor_rec(d, factors);
    factor_rec(n / d, factors);
}

u64 largest_prime_factor(u64 n) {
    if (n < 2) return 1;
    if (is_prime(n)) return n;

    std::vector<u64> factors;
    factors.reserve(8);
    factor_rec(n, factors);
    return *std::max_element(factors.begin(), factors.end());
}

std::vector<int> sieve_primes(int limit) {
    std::vector<char> is(limit + 1, true);
    is[0] = false;
    is[1] = false;
    for (int i = 2; 1LL * i * i <= limit; ++i) {
        if (!is[i]) continue;
        for (int j = i * i; j <= limit; j += i) is[j] = false;
    }

    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (is[i]) primes.push_back(i);
    }
    return primes;
}

u64 compute_h_small(u64 n) {
    if (n < 2) return 0;

    const u64 block_size = 5000000ULL;
    const int sq = static_cast<int>(std::sqrt(static_cast<long double>(n + 8))) + 5;
    const std::vector<int> primes = sieve_primes(sq);

    u64 best = 0;
    std::vector<u64> tail;
    tail.reserve(8);

    for (u64 l = 2; l <= n + 8; l += block_size) {
        const u64 r = std::min(n + 8, l + block_size - 1);
        const size_t len = static_cast<size_t>(r - l + 1);

        std::vector<u64> rem(len);
        std::vector<u64> lpf(len, 1);
        for (size_t i = 0; i < len; ++i) rem[i] = l + static_cast<u64>(i);

        for (int p : primes) {
            const u64 pp = static_cast<u64>(p);
            if (pp * pp > r) break;
            u64 start = ((l + pp - 1) / pp) * pp;
            for (u64 m = start; m <= r; m += pp) {
                const size_t idx = static_cast<size_t>(m - l);
                if (rem[idx] % pp != 0) continue;
                lpf[idx] = pp;
                do {
                    rem[idx] /= pp;
                } while (rem[idx] % pp == 0);
            }
        }

        for (size_t i = 0; i < len; ++i) {
            if (rem[i] > 1) lpf[i] = std::max(lpf[i], rem[i]);
        }

        std::vector<u64> merged;
        merged.reserve(tail.size() + len);
        merged.insert(merged.end(), tail.begin(), tail.end());
        merged.insert(merged.end(), lpf.begin(), lpf.end());

        const u64 base_n = l - static_cast<u64>(tail.size());
        if (merged.size() >= 9) {
            u64 window = 0;
            for (int i = 0; i < 9; ++i) window += merged[i];
            for (size_t j = 0; j + 8 < merged.size(); ++j) {
                const u64 cur_n = base_n + static_cast<u64>(j);
                if (cur_n >= 2 && cur_n <= n && window > best) best = window;
                if (j + 9 < merged.size()) {
                    window += merged[j + 9];
                    window -= merged[j];
                }
            }
        }

        tail.clear();
        const size_t take = std::min<size_t>(8, merged.size());
        tail.insert(tail.end(), merged.end() - take, merged.end());
    }

    return best;
}

struct Best {
    std::atomic<u64> g;
    std::atomic<u64> n;
};

void update_best(Best& best, u64 n, u64 g) {
    u64 cur = best.g.load(std::memory_order_relaxed);
    while (g > cur && !best.g.compare_exchange_weak(cur, g, std::memory_order_relaxed)) {}
    if (g == best.g.load(std::memory_order_relaxed)) {
        best.n.store(n, std::memory_order_relaxed);
    }
}

bool validate_examples() {
    if (largest_prime_factor(100) != 5ULL) {
        std::cerr << "Validation failed: f(100) != 5\n";
        return false;
    }
    if (largest_prime_factor(101) != 101ULL) {
        std::cerr << "Validation failed: f(101) != 101\n";
        return false;
    }

    u64 g100 = 0;
    for (u64 x = 100; x <= 108; ++x) g100 += largest_prime_factor(x);
    if (g100 != 409ULL) {
        std::cerr << "Validation failed: g(100) != 409\n";
        return false;
    }

    if (compute_h_small(100ULL) != 417ULL) {
        std::cerr << "Validation failed: h(100) != 417\n";
        return false;
    }

    const u64 h1e9 = compute_h_small(1000000000ULL);
    if (h1e9 != 4896292593ULL) {
        std::cerr << "Validation failed: h(1e9) mismatch, got " << h1e9 << "\n";
        return false;
    }

    return true;
}

struct ResidueInfo {
    int r;
    std::array<int, 9> d;
    long double K;
};

u64 solve_h_1e16() {
    Best best{kSeedBestG, kSeedBestN};

    // Sanity-check the seed lower bound.
    {
        const u64 n = kSeedBestN;
        u64 g = 0;
        for (int i = 0; i < 9; ++i) g += largest_prime_factor(n + static_cast<u64>(i));
        if (g != kSeedBestG) {
            std::cerr << "Internal error: seed bound mismatch\n";
            std::exit(1);
        }
    }

    // Prove all non-top residues are already below the seed bound.
    for (int r = 0; r < 2520; ++r) {
        if (r == 311 || r == 2201) continue;
        long double C = 0;
        long double K = 0;
        for (int i = 0; i < 9; ++i) {
            const int d = std::gcd(2520, r + i);
            C += 1.0L / static_cast<long double>(d);
            K += static_cast<long double>(i) / static_cast<long double>(d);
        }
        const long double ub = C * static_cast<long double>(kTargetN) + K;
        if (ub > static_cast<long double>(kSeedBestG)) {
            std::cerr << "Internal error: residue " << r
                      << " not eliminated by seed lower bound\n";
            std::exit(1);
        }
    }

    const std::array<ResidueInfo, 2> infos = {
        ResidueInfo{311,  {1, 24, 1, 2, 315, 4, 1, 6, 1}, kK311},
        ResidueInfo{2201, {1, 6, 1, 4, 315, 2, 1, 24, 1}, kK2201},
    };

    // Candidate n must satisfy n > best/C_top.
    const u64 n_start = static_cast<u64>(static_cast<long double>(kSeedBestG) / kCTop) + 1;

    // Prime-filter safety: dropping any mandatory prime term loses more than
    // the maximal remaining UB margin in [n_start, N].
    for (const auto& info : infos) {
        const long double diff = (kCTop * static_cast<long double>(kTargetN) + info.K) -
                                 static_cast<long double>(kSeedBestG);

        auto loss = [&](int i, int m) {
            const long double d = static_cast<long double>(info.d[i]);
            const long double coeff = (1.0L / d) - (1.0L / (d * static_cast<long double>(m)));
            const long double constant = (static_cast<long double>(i) / d) -
                                         (static_cast<long double>(i) / (d * static_cast<long double>(m)));
            return coeff * static_cast<long double>(n_start) + constant;
        };

        // Odd terms (d=1), composite implies a prime factor >=11.
        for (int i : {0, 2, 6, 8}) {
            if (loss(i, 11) <= diff) {
                std::cerr << "Internal error: odd-prime filter not safe\n";
                std::exit(1);
            }
        }

        if (info.r == 311) {
            if (loss(1, 11) <= diff || loss(3, 11) <= diff ||
                loss(5, 11) <= diff || loss(7, 11) <= diff) {
                std::cerr << "Internal error: quotient filter not safe for residue 311\n";
                std::exit(1);
            }
        } else {
            if (loss(1, 11) <= diff || loss(3, 11) <= diff || loss(5, 11) <= diff ||
                loss(7, 2) <= diff) {
                std::cerr << "Internal error: quotient filter not safe for residue 2201\n";
                std::exit(1);
            }
        }
    }

    std::array<std::atomic<u64>, 2> next_t;
    std::array<u64, 2> t_min{};

    for (int idx = 0; idx < 2; ++idx) {
        const u64 r = static_cast<u64>(infos[idx].r);
        const u64 t_max = (kTargetN - r) / kStep;
        next_t[idx].store(t_max, std::memory_order_relaxed);
        t_min[idx] = (n_start <= r) ? 0 : ((n_start - r + kStep - 1) / kStep);
    }

    const u64 kChunk = 200000ULL;
    std::atomic<u64> checked{0};
    std::atomic<u64> odd_hits{0};
    std::atomic<u64> hit8{0};
    std::mutex print_mu;

    auto worker = [&](int worker_id) {
        const int idx = worker_id & 1;
        const auto& info = infos[idx];

        while (true) {
            const u64 hi = next_t[idx].fetch_sub(kChunk, std::memory_order_relaxed);
            if (hi < t_min[idx]) break;

            u64 lo = (hi >= kChunk) ? (hi - kChunk + 1) : 0;
            if (lo < t_min[idx]) lo = t_min[idx];
            if (lo > hi) continue;

            for (u64 t = hi;; --t) {
                const u64 n = kStep * t + static_cast<u64>(info.r);

                const long double ub = kCTop * static_cast<long double>(n) + info.K;
                if (ub <= static_cast<long double>(best.g.load(std::memory_order_relaxed))) break;

                checked.fetch_add(1, std::memory_order_relaxed);

                if (!is_prime(n) || !is_prime(n + 2) || !is_prime(n + 6) || !is_prime(n + 8)) {
                    if (t == lo) break;
                    continue;
                }
                odd_hits.fetch_add(1, std::memory_order_relaxed);

                u64 q1 = 0, q3 = 0, q5 = 0, q7 = 0;
                if (info.r == 311) {
                    q1 = (n + 1) / 24;
                    q3 = (n + 3) / 2;
                    q5 = (n + 5) / 4;
                    q7 = (n + 7) / 6;
                } else {
                    q1 = (n + 1) / 6;
                    q3 = (n + 3) / 4;
                    q5 = (n + 5) / 2;
                    q7 = (n + 7) / 24;
                }

                if (!is_prime(q1) || !is_prime(q3) || !is_prime(q5) || !is_prime(q7)) {
                    if (t == lo) break;
                    continue;
                }
                hit8.fetch_add(1, std::memory_order_relaxed);

                const u64 term4 = largest_prime_factor(n + 4);
                const u64 g = n + (n + 2) + (n + 6) + (n + 8) + q1 + q3 + q5 + q7 + term4;

                if (g > best.g.load(std::memory_order_relaxed)) {
                    update_best(best, n, g);
                    std::lock_guard<std::mutex> lk(print_mu);
                    std::cout << "new best n=" << n << " g=" << g
                              << " ratio=" << std::setprecision(12)
                              << (static_cast<long double>(g) / static_cast<long double>(n))
                              << " residue=" << info.r << '\n';
                }

                if (t == lo) break;
            }
        }
    };

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 8;
    threads = std::min<unsigned>(threads, 8);

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned i = 0; i < threads; ++i) pool.emplace_back(worker, static_cast<int>(i));
    for (auto& th : pool) th.join();

    std::cout << "Large-scan stats: checked=" << checked.load()
              << " odd_quadruplets=" << odd_hits.load()
              << " full_8prime_hits=" << hit8.load() << '\n';

    return best.g.load();
}

}  // namespace

int main() {
    if (!validate_examples()) return 1;

    const u64 ans = solve_h_1e16();
    std::cout << ans << '\n';
    return 0;
}
