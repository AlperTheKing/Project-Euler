#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using int64 = long long;
using i128 = __int128_t;
using u128 = __uint128_t;

static constexpr int64 kMod = 1000000000LL;
static constexpr int64 kDefaultN = 100000000000LL;
static constexpr int kSieveLimit = 20000000;

static inline int64 mod_norm(int64 x) {
    x %= kMod;
    if (x < 0) x += kMod;
    return x;
}

static inline int64 mod_mul(int64 a, int64 b) {
    return static_cast<int64>((static_cast<i128>(a) * b) % kMod);
}

static inline int64 sum_arith_mod(int64 l, int64 r) {
    u128 t = static_cast<u128>(l + r) * static_cast<u128>(r - l + 1) / 2;
    return static_cast<int64>(t % kMod);
}

// A(M) = sum_{n<=M} sigma(n) = sum_{k<=M} k * floor(M/k).
static int64 sum_sigma_prefix(int64 M) {
    int64 total = 0;
    for (int64 l = 1, r; l <= M; l = r + 1) {
        int64 q = M / l;
        r = M / q;
        u128 t = static_cast<u128>(q) * static_cast<u128>(q + 1) / 2;
        t *= static_cast<u128>(r - l + 1);
        total += static_cast<int64>(t % kMod);
        if (total >= kMod) total %= kMod;
    }
    return total % kMod;
}

// Summatory of k * mu(k) with memoization.
class MobiusPrefix {
public:
    explicit MobiusPrefix(int64 max_n, size_t cache_hint) {
        limit_ = static_cast<int>(max_n);
        mu_.assign(limit_ + 1, 0);
        prefix_.assign(limit_ + 1, 0);
        cache_.reserve(cache_hint);
        build();
    }

    int64 sum_prefix(int64 n) {
        if (n <= 0) return 0;
        if (n <= limit_) return prefix_[static_cast<size_t>(n)];
        auto it = cache_.find(n);
        if (it != cache_.end()) return it->second;

        int64 ans = 1;
        for (int64 l = 2, r; l <= n; l = r + 1) {
            int64 v = n / l;
            r = n / v;
            int64 sum_lr = sum_arith_mod(l, r);
            int64 sub = mod_mul(sum_lr, sum_prefix(v));
            ans -= sub;
            if (ans < 0) ans += kMod;
        }
        cache_[n] = static_cast<int32_t>(ans);
        return ans;
    }

private:
    int limit_ = 0;
    std::vector<int8_t> mu_;
    std::vector<int32_t> prefix_;
    std::unordered_map<int64, int32_t> cache_;

    void build() {
        std::vector<int> primes;
        primes.reserve(limit_ / 10);
        std::vector<uint8_t> is_comp(limit_ + 1, 0);
        mu_[1] = 1;
        for (int i = 2; i <= limit_; ++i) {
            if (!is_comp[i]) {
                primes.push_back(i);
                mu_[i] = -1;
            }
            for (int p : primes) {
                int64 v = static_cast<int64>(i) * p;
                if (v > limit_) break;
                is_comp[static_cast<size_t>(v)] = 1;
                if (i % p == 0) {
                    mu_[static_cast<size_t>(v)] = 0;
                    break;
                }
                mu_[static_cast<size_t>(v)] = static_cast<int8_t>(-mu_[i]);
            }
        }
        int64 running = 0;
        for (int i = 1; i <= limit_; ++i) {
            running += static_cast<int64>(i) * mu_[i];
            running = mod_norm(running);
            prefix_[static_cast<size_t>(i)] = static_cast<int32_t>(running);
        }
    }
};

struct RangeData {
    int64 l;
    int64 r;
    int64 m;
    int64 mu_sum;
};

// Uses identity:
//   sigma(xy) = sum_{d|gcd(x,y)} mu(d) * d * sigma(x/d) * sigma(y/d).
// Summing over i,j<=N gives S(N)=sum_{d<=N} mu(d)*d*A(N/d)^2.
static int64 solve(int64 N, MobiusPrefix& pref, int threads) {
    std::vector<RangeData> ranges;
    ranges.reserve(static_cast<size_t>(2.0 * std::sqrt(static_cast<long double>(N)) + 10));
    for (int64 l = 1, r; l <= N; l = r + 1) {
        int64 q = N / l;
        r = N / q;
        ranges.push_back({l, r, q, 0});
    }

    int64 prev_prefix = 0;
    for (auto& range : ranges) {
        int64 pref_r = pref.sum_prefix(range.r);
        int64 pref_lm1 = (range.l == 1) ? 0 : prev_prefix;
        int64 mu_sum = pref_r - pref_lm1;
        if (mu_sum < 0) mu_sum += kMod;
        range.mu_sum = mu_sum;
        prev_prefix = pref_r;
    }

    if (threads < 1) threads = 1;
    threads = std::min<int>(threads, static_cast<int>(ranges.size()));
    std::atomic<size_t> next_idx(0);
    std::vector<int64> partial(static_cast<size_t>(threads), 0);

    auto worker = [&](int tid) {
        int64 local = 0;
        while (true) {
            size_t idx = next_idx.fetch_add(1, std::memory_order_relaxed);
            if (idx >= ranges.size()) break;
            const auto& range = ranges[idx];
            if (range.mu_sum == 0) continue;
            int64 a = sum_sigma_prefix(range.m);
            int64 term = mod_mul(range.mu_sum, mod_mul(a, a));
            local += term;
            if (local >= kMod) local %= kMod;
        }
        partial[static_cast<size_t>(tid)] = local % kMod;
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(threads));
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto& th : pool) th.join();

    int64 total = 0;
    for (int64 v : partial) {
        total += v;
        if (total >= kMod) total %= kMod;
    }
    return total % kMod;
}

static bool validate(MobiusPrefix& pref, int threads) {
    struct Test {
        int64 N;
        int64 expected_mod;
    };
    const std::vector<Test> tests = {
        {3, 59},
        {1000, 563576517282LL % kMod},
        {100000, 215766508},
    };
    bool ok = true;
    for (const auto& test : tests) {
        int64 got = solve(test.N, pref, threads);
        std::cout << "Validation S(" << test.N << ") mod 1e9: " << got;
        if (got == test.expected_mod) {
            std::cout << " [PASS]";
        } else {
            std::cout << " [FAIL] Expected " << test.expected_mod;
            ok = false;
        }
        std::cout << "\n";
    }
    return ok;
}

int main(int argc, char** argv) {
    int64 N = kDefaultN;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads == 0) threads = 4;
    bool run_validation = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--threads=", 0) == 0) {
            threads = std::max(1, std::stoi(arg.substr(10)));
        } else if (arg == "--no-validate") {
            run_validation = false;
        } else {
            N = std::stoll(arg);
        }
    }

    int64 sieve_limit = std::min<int64>(kSieveLimit, N);
    size_t cache_hint = static_cast<size_t>(2.0 * std::sqrt(static_cast<long double>(N)) + 1000);
    MobiusPrefix pref(sieve_limit, cache_hint);

    if (run_validation) {
        std::cout << "Running validations...\n";
        if (!validate(pref, std::min(threads, 4))) {
            std::cout << "Validation failed. Aborting.\n";
            return 1;
        }
        std::cout << "--------------------------------\n";
    }

    std::cout << "Solving for N = " << N << "...\n";
    auto start = std::chrono::high_resolution_clock::now();
    int64 result = solve(N, pref, threads);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Result S(N) mod 1e9: " << result << "\n";
    std::cout << "Time elapsed: " << elapsed.count() << "s\n";
    return 0;
}
