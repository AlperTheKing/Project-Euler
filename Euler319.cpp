#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

using int64 = long long;
using i128 = __int128_t;

namespace {
constexpr int64 kMod = 1000000000LL;
constexpr int64 kMod2 = 2 * kMod;

int64 mod_norm(int64 x) {
    x %= kMod;
    if (x < 0) x += kMod;
    return x;
}

int64 mul_mod(int64 a, int64 b, int64 mod) {
    return static_cast<int64>((static_cast<__int128>(a) * b) % mod);
}

int64 pow_mod(int64 base, int64 exp, int64 mod) {
    int64 res = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = mul_mod(res, base, mod);
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }
    return res;
}

int64 sum_a_prefix(int64 m) {
    if (m <= 0) return 0;
    int64 pow3 = pow_mod(3, m + 1, kMod2);
    int64 num = pow3 - 3;
    if (num < 0) num += kMod2;
    int64 sum3 = num / 2;  // (3^{m+1}-3)/2 mod 1e9

    int64 pow2 = pow_mod(2, m + 1, kMod);
    int64 sum2 = pow2 - 2;
    if (sum2 < 0) sum2 += kMod;

    int64 res = sum3 - sum2 - (m % kMod);
    return mod_norm(res);
}

class Mertens {
public:
    explicit Mertens(int64 n) {
        limit_ = static_cast<int64>(std::pow(static_cast<long double>(n), 2.0L / 3.0L)) + 1;
        if (limit_ < 1) limit_ = 1;
        init_mu();
        cache_.reserve(1 << 20);
    }

    int64 get(int64 n) {
        if (n <= limit_) return prefix_[static_cast<size_t>(n)];
        auto it = cache_.find(n);
        if (it != cache_.end()) return it->second;
        int64 res = 1;
        int64 l = 2;
        while (l <= n) {
            int64 q = n / l;
            int64 r = n / q;
            res -= (r - l + 1) * get(q);
            l = r + 1;
        }
        cache_[n] = res;
        return res;
    }

private:
    int64 limit_ = 0;
    std::vector<int> mu_;
    std::vector<int> prefix_;
    std::unordered_map<int64, int64> cache_;

    void init_mu() {
        mu_.assign(static_cast<size_t>(limit_) + 1, 0);
        prefix_.assign(static_cast<size_t>(limit_) + 1, 0);
        std::vector<int> primes;
        std::vector<int> is_comp(static_cast<size_t>(limit_) + 1, 0);
        mu_[1] = 1;
        for (int i = 2; i <= limit_; ++i) {
            if (!is_comp[static_cast<size_t>(i)]) {
                primes.push_back(i);
                mu_[static_cast<size_t>(i)] = -1;
            }
            for (int p : primes) {
                int64 v = static_cast<int64>(i) * p;
                if (v > limit_) break;
                is_comp[static_cast<size_t>(v)] = 1;
                if (i % p == 0) {
                    mu_[static_cast<size_t>(v)] = 0;
                    break;
                } else {
                    mu_[static_cast<size_t>(v)] = -mu_[static_cast<size_t>(i)];
                }
            }
        }
        for (int i = 1; i <= limit_; ++i) {
            prefix_[static_cast<size_t>(i)] = prefix_[static_cast<size_t>(i - 1)] + mu_[static_cast<size_t>(i)];
        }
    }
};

struct Segment {
    int64 l;
    int64 r;
    int64 q;
    int64 mq;
};

int64 compute_t_mod(int64 n, int threads) {
    Mertens mertens(n);
    std::vector<Segment> segs;
    segs.reserve(static_cast<size_t>(2 * std::sqrt(static_cast<long double>(n)) + 10));

    for (int64 l = 1; l <= n; ) {
        int64 q = n / l;
        int64 r = n / q;
        segs.push_back({l, r, q, 0});
        l = r + 1;
    }

    for (auto& seg : segs) {
        seg.mq = mertens.get(seg.q);
    }

    if (threads <= 0) threads = 1;
    threads = std::min<int>(threads, static_cast<int>(segs.size()));

    auto worker = [&](size_t start, size_t step) -> int64 {
        int64 local = 0;
        for (size_t i = start; i < segs.size(); i += step) {
            const auto& seg = segs[i];
            int64 sum_a = sum_a_prefix(seg.r) - sum_a_prefix(seg.l - 1);
            sum_a = mod_norm(sum_a);
            int64 mq_mod = seg.mq % kMod;
            if (mq_mod < 0) mq_mod += kMod;
            int64 contrib = static_cast<int64>((static_cast<__int128>(sum_a) * mq_mod) % kMod);
            local += contrib;
            if (local >= kMod || local <= -kMod) local %= kMod;
        }
        return mod_norm(local);
    };

    int64 total = 1 % kMod;
    if (threads == 1 || segs.size() < 1024) {
        total = mod_norm(total + worker(0, 1));
        return total;
    }

    std::vector<int64> partial(static_cast<size_t>(threads), 0);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(threads));
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            partial[static_cast<size_t>(t)] = worker(static_cast<size_t>(t), static_cast<size_t>(threads));
        });
    }
    for (auto& th : pool) th.join();

    for (int64 val : partial) {
        total += val;
        if (total >= kMod || total <= -kMod) total %= kMod;
    }
    return mod_norm(total);
}

int64 pow_ll(int64 base, int64 exp) {
    int64 res = 1;
    while (exp > 0) {
        if (exp & 1) res *= base;
        base *= base;
        exp >>= 1;
    }
    return res;
}

int64 t_direct(int n) {
    std::vector<int> mu(n + 1, 0);
    std::vector<int> is_comp(n + 1, 0);
    std::vector<int> primes;
    mu[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (!is_comp[i]) {
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            int v = i * p;
            if (v > n) break;
            is_comp[v] = 1;
            if (i % p == 0) {
                mu[v] = 0;
                break;
            } else {
                mu[v] = -mu[i];
            }
        }
    }

    std::vector<int> pref(n + 1, 0);
    for (int i = 1; i <= n; ++i) pref[i] = pref[i - 1] + mu[i];

    i128 total = 1;
    for (int k = 1; k <= n; ++k) {
        i128 a = static_cast<i128>(pow_ll(3, k)) - pow_ll(2, k) - 1;
        total += a * pref[n / k];
    }
    return static_cast<int64>(total);
}

void run_validation() {
    struct Test {
        int n;
        int64 expected;
    } tests[] = {
        {2, 5},
        {5, 293},
        {10, 86195},
        {20, 5227991891LL},
    };

    for (const auto& test : tests) {
        int64 val = t_direct(test.n);
        if (val != test.expected) {
            std::cerr << "Validation failed for n=" << test.n << ": got " << val
                      << ", expected " << test.expected << "\n";
            std::exit(1);
        }
        int64 mod_val = compute_t_mod(test.n, 1);
        int64 exp_mod = test.expected % kMod;
        if (mod_val != exp_mod) {
            std::cerr << "Mod validation failed for n=" << test.n << ": got " << mod_val
                      << ", expected " << exp_mod << "\n";
            std::exit(1);
        }
    }
}

}  // namespace

int main() {
    run_validation();

    const int64 n = 10000000000LL;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads <= 0) threads = 1;
    int64 ans = compute_t_mod(n, threads);
    std::cout << ans << "\n";
    return 0;
}
