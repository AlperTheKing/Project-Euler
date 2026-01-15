#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using int64 = long long;
using u128 = unsigned __int128;

constexpr int64 MOD = 1000000007LL;

struct Range {
    int64 q = 0;
    u128 count = 0;
    int64 count_mod = 0;
};

struct Partial {
    int64 h = 1;
    int64 s1 = 0;
    int64 s2 = 0;
};

int64 mul_mod(int64 a, int64 b) {
    return static_cast<int64>((__int128)a * b % MOD);
}

int64 mod_pow(int64 base, u128 exp) {
    if (exp == 0) return 1;
    base %= MOD;
    if (base == 0) return 0;
    int64 result = 1;
    while (exp > 0) {
        if (exp & 1) result = mul_mod(result, base);
        base = mul_mod(base, base);
        exp >>= 1;
    }
    return result;
}

int64 mod_u128(u128 value) {
    return static_cast<int64>(value % MOD);
}

class PhiSummatory {
  public:
    explicit PhiSummatory(int64 n) { init(n); }

    u128 sum(int64 n) {
        if (n <= 0) return 0;
        if (n <= limit_) return prefix_[static_cast<std::size_t>(n)];
        auto it = memo_.find(n);
        if (it != memo_.end()) return it->second;
        u128 res = static_cast<u128>(n) * static_cast<u128>(n + 1) / 2;
        int64 i = 2;
        while (i <= n) {
            int64 q = n / i;
            int64 j = n / q;
            res -= static_cast<u128>(j - i + 1) * sum(q);
            i = j + 1;
        }
        memo_[n] = res;
        return res;
    }

  private:
    int64 limit_ = 0;
    std::vector<uint64_t> prefix_;
    std::unordered_map<int64, u128> memo_;

    void init(int64 n) {
        int64 limit = static_cast<int64>(std::pow(static_cast<long double>(n), 2.0L / 3.0L));
        const int64 kMinLimit = 1000000;
        const int64 kMaxLimit = 10000000;
        if (limit < kMinLimit) limit = std::min<int64>(n, kMinLimit);
        if (limit > kMaxLimit) limit = kMaxLimit;
        if (limit < 1) limit = 1;
        if (limit > n) limit = n;
        limit_ = limit;

        std::vector<uint32_t> phi(limit_ + 1, 0);
        std::vector<int> primes;
        primes.reserve(static_cast<std::size_t>(limit_ / 10 + 10));
        if (limit_ >= 1) phi[1] = 1;
        for (int i = 2; i <= limit_; ++i) {
            if (phi[i] == 0) {
                phi[i] = static_cast<uint32_t>(i - 1);
                primes.push_back(i);
            }
            for (int p : primes) {
                long long v = 1LL * p * i;
                if (v > limit_) break;
                if (i % p == 0) {
                    phi[static_cast<std::size_t>(v)] = phi[i] * static_cast<uint32_t>(p);
                    break;
                }
                phi[static_cast<std::size_t>(v)] = phi[i] * static_cast<uint32_t>(p - 1);
            }
        }

        prefix_.assign(limit_ + 1, 0);
        for (int64 i = 1; i <= limit_; ++i) {
            prefix_[static_cast<std::size_t>(i)] =
                prefix_[static_cast<std::size_t>(i - 1)] + phi[static_cast<std::size_t>(i)];
        }
        std::vector<uint32_t>().swap(phi);
        memo_.reserve(1 << 20);
    }
};

int64 default_n() {
    int64 n = 1;
    for (int i = 0; i < 13; ++i) n *= 7;
    return n;
}

std::vector<Range> build_ranges(int64 n, PhiSummatory& phi) {
    std::vector<Range> ranges;
    int64 approx = static_cast<int64>(2 * std::sqrt(static_cast<long double>(n)) + 10);
    ranges.reserve(static_cast<std::size_t>(approx));
    for (int64 r = 1; r <= n;) {
        int64 q = n / r;
        int64 R = n / q;
        u128 sum_phi = phi.sum(R) - phi.sum(r - 1);
        u128 count = sum_phi * 4;
        ranges.push_back(Range{q, count, mod_u128(count)});
        r = R + 1;
    }
    return ranges;
}

Partial compute_partial(const std::vector<Range>& ranges, std::size_t start, std::size_t end) {
    Partial out;
    for (std::size_t i = start; i < end; ++i) {
        const Range& range = ranges[i];
        int64 base = (range.q + 1) % MOD;
        out.h = mul_mod(out.h, mod_pow(base, range.count));
        int64 q_mod = range.q % MOD;
        int64 q2 = mul_mod(q_mod, q_mod);
        out.s1 = (out.s1 + mul_mod(range.count_mod, q_mod)) % MOD;
        out.s2 = (out.s2 + mul_mod(range.count_mod, q2)) % MOD;
    }
    return out;
}

// Directions come in opposite pairs with the same weight w.
// Total weighted subsets are H^2 with H = product over pairs (1 + w).
// Bad subsets in a closed half-plane contribute B = 1 + 2*H*sum(w) - sum(w^2),
// so P = H^2 - B.
int64 compute_polar_polygons(int64 n, unsigned threads) {
    if (n <= 0) return 0;
    PhiSummatory phi(n);
    auto ranges = build_ranges(n, phi);

    if (threads == 0) threads = 1;
    unsigned use_threads = std::min<unsigned>(threads, ranges.size());
    if (use_threads <= 1) {
        auto partial = compute_partial(ranges, 0, ranges.size());
        int64 h = partial.h;
        int64 s1 = partial.s1;
        int64 s2 = partial.s2;
        int64 result = (mul_mod(h, h) - 1 - mul_mod((2 * h) % MOD, s1) + s2) % MOD;
        if (result < 0) result += MOD;
        return result;
    }

    std::vector<std::thread> pool;
    std::vector<Partial> partials(use_threads);
    std::size_t chunk = (ranges.size() + use_threads - 1) / use_threads;
    pool.reserve(use_threads);

    for (unsigned t = 0; t < use_threads; ++t) {
        std::size_t start = t * chunk;
        std::size_t end = std::min(ranges.size(), start + chunk);
        if (start >= end) continue;
        pool.emplace_back([&, start, end, t]() {
            partials[t] = compute_partial(ranges, start, end);
        });
    }
    for (auto& th : pool) th.join();

    int64 h = 1;
    int64 s1 = 0;
    int64 s2 = 0;
    for (const auto& part : partials) {
        h = mul_mod(h, part.h);
        s1 += part.s1;
        if (s1 >= MOD) s1 -= MOD;
        s2 += part.s2;
        if (s2 >= MOD) s2 -= MOD;
    }

    int64 result = (mul_mod(h, h) - 1 - mul_mod((2 * h) % MOD, s1) + s2) % MOD;
    if (result < 0) result += MOD;
    return result;
}

bool validate() {
    if (compute_polar_polygons(1, 1) != 131) {
        std::cerr << "Validation failed: P(1) mismatch.\n";
        return false;
    }
    if (compute_polar_polygons(2, 1) != 1648531) {
        std::cerr << "Validation failed: P(2) mismatch.\n";
        return false;
    }
    if (compute_polar_polygons(3, 1) != 461288482) {
        std::cerr << "Validation failed: P(3) mismatch (mod).\n";
        return false;
    }
    if (compute_polar_polygons(343, 1) != 937293740) {
        std::cerr << "Validation failed: P(343) mismatch (mod).\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    int64 n = default_n();
    if (argc >= 2) {
        n = std::stoll(argv[1]);
    }
    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0) threads = 1;
    if (argc >= 3) {
        threads = static_cast<unsigned>(std::stoul(argv[2]));
        if (threads == 0) threads = 1;
    }

    if (!validate()) return 1;

    int64 result = compute_polar_polygons(n, threads);
    std::cout << result << '\n';
    return 0;
}
