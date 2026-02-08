#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = __uint128_t;

constexpr u64 kMod = 1'000'000'007ULL;
constexpr u64 kDefaultN = 1'000'000'000'000ULL;
constexpr u64 kCheckpointN1 = 10ULL;
constexpr u64 kCheckpointExpected1 = 3'053ULL;
constexpr u64 kCheckpointN2 = 100'000ULL;
constexpr u64 kCheckpointExpected2 = 157'612'967ULL;

struct Options {
    u64 n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct U64Hash {
    std::size_t operator()(const u64 value) const noexcept {
        u64 x = value + 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
        x ^= (x >> 31U);
        return static_cast<std::size_t>(x);
    }
};

bool parse_u64_after_prefix(const std::string& arg, const std::string& prefix, u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }
    const std::string tail = arg.substr(prefix.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(ch - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const std::string& prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--n=", options.n)) {
            continue;
        }
        if (parse_unsigned_after_prefix(arg, "--threads=", options.requested_threads)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

unsigned pick_thread_count(const bool allow_multithreading, const unsigned requested_threads) {
    if (!allow_multithreading) {
        return 1U;
    }
    if (requested_threads > 0U) {
        return std::max(1U, requested_threads);
    }
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0U) {
        hw = 4U;
    }
    return hw;
}

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

inline u64 add_mod(const u64 a, const u64 b) {
    const u64 s = a + b;
    return (s >= kMod) ? (s - kMod) : s;
}

inline u64 sub_mod(const u64 a, const u64 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

inline u64 mul_mod(const u64 a, const u64 b) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % static_cast<u128>(kMod));
}

int chi4(const u64 n) {
    if ((n & 1ULL) == 0ULL) {
        return 0;
    }
    return ((n & 3ULL) == 1ULL) ? 1 : -1;
}

u64 sum_cubes_1_to_mod(const u64 n) {
    const u128 a = (static_cast<u128>(n) * static_cast<u128>(n + 1ULL)) / static_cast<u128>(2ULL);
    const u64 am = static_cast<u64>(a % static_cast<u128>(kMod));
    return mul_mod(am, am);
}

u64 prefix_chi_1_to_mod(const u64 n) {
    if (n == 0ULL) {
        return 0ULL;
    }
    const u64 r = n & 3ULL;
    return ((r == 1ULL) || (r == 2ULL)) ? 1ULL : 0ULL;
}

std::vector<int> sieve_primes(const int limit) {
    if (limit < 2) {
        return {};
    }

    const int size = limit / 2 + 1;
    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(size), 1U);
    is_prime[0] = 0U;  // 1 is not prime

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(limit)));
    for (int i = 1; (2 * i + 1) <= root; ++i) {
        if (is_prime[static_cast<std::size_t>(i)] == 0U) {
            continue;
        }
        const int p = 2 * i + 1;
        int start = (p * p) / 2;
        for (int j = start; j < size; j += p) {
            is_prime[static_cast<std::size_t>(j)] = 0U;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::max(1.0L, std::log(static_cast<long double>(limit)))));
    primes.push_back(2);
    for (int i = 1; i < size; ++i) {
        if (is_prime[static_cast<std::size_t>(i)] != 0U) {
            primes.push_back(2 * i + 1);
        }
    }
    return primes;
}

struct ValuesIndex {
    std::vector<u64> values;
    u64 root = 0ULL;
    std::vector<int> idx_small;
    std::vector<int> idx_large;

    int index_of(const u64 n, const u64 x) const {
        if (x <= root) {
            return idx_small[static_cast<std::size_t>(x)];
        }
        return idx_large[static_cast<std::size_t>(n / x)];
    }
};

ValuesIndex build_values_and_index(const u64 n) {
    ValuesIndex out;
    out.root = isqrt_u64(n);

    std::vector<u64> large;
    large.reserve(static_cast<std::size_t>(2ULL * out.root + 16ULL));
    for (u64 i = 1ULL; i <= n;) {
        const u64 v = n / i;
        large.push_back(v);
        i = n / v + 1ULL;
    }

    std::vector<u64> small;
    small.reserve(static_cast<std::size_t>(out.root));
    for (u64 v = out.root; v >= 1ULL; --v) {
        small.push_back(v);
    }

    out.values.reserve(large.size() + small.size());
    std::size_t ia = 0U;
    std::size_t ib = 0U;
    while (ia < large.size() && ib < small.size()) {
        const u64 a = large[ia];
        const u64 b = small[ib];
        if (a > b) {
            out.values.push_back(a);
            ++ia;
        } else if (a < b) {
            out.values.push_back(b);
            ++ib;
        } else {
            out.values.push_back(a);
            ++ia;
            ++ib;
        }
    }
    while (ia < large.size()) {
        out.values.push_back(large[ia++]);
    }
    while (ib < small.size()) {
        out.values.push_back(small[ib++]);
    }

    out.idx_small.assign(static_cast<std::size_t>(out.root + 1ULL), 0);
    out.idx_large.assign(static_cast<std::size_t>(out.root + 1ULL), 0);

    for (std::size_t idx = 0; idx < out.values.size(); ++idx) {
        const u64 v = out.values[idx];
        if (v <= out.root) {
            out.idx_small[static_cast<std::size_t>(v)] = static_cast<int>(idx);
        } else {
            out.idx_large[static_cast<std::size_t>(n / v)] = static_cast<int>(idx);
        }
    }

    return out;
}

void parallel_chunked_for(const unsigned threads,
                          const std::size_t begin,
                          const std::size_t end,
                          const std::function<void(std::size_t, std::size_t)>& work) {
    if (threads <= 1U || end <= begin + 1U) {
        work(begin, end);
        return;
    }

    const std::size_t length = end - begin;
    const unsigned use_threads = std::min<unsigned>(threads, static_cast<unsigned>(length));
    if (use_threads <= 1U) {
        work(begin, end);
        return;
    }

    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(use_threads));

    std::size_t chunk_begin = begin;
    for (unsigned t = 0U; t < use_threads; ++t) {
        const std::size_t remaining = end - chunk_begin;
        const std::size_t lanes = static_cast<std::size_t>(use_threads - t);
        const std::size_t chunk_size = (remaining + lanes - 1U) / lanes;
        const std::size_t chunk_end = chunk_begin + chunk_size;
        pool.emplace_back([&, chunk_begin, chunk_end]() { work(chunk_begin, chunk_end); });
        chunk_begin = chunk_end;
    }

    for (std::thread& th : pool) {
        th.join();
    }
}

struct PrimeSumData {
    ValuesIndex index;
    std::vector<u64> gprime;
};

PrimeSumData compute_prime_sums_gprime(const u64 n,
                                       const std::vector<int>& primes,
                                       const unsigned threads) {
    PrimeSumData out;
    out.index = build_values_and_index(n);
    const std::size_t m = out.index.values.size();

    std::vector<u64> gcube(m, 0ULL);
    std::vector<u64> gchi(m, 0ULL);

    parallel_chunked_for(threads, 0U, m, [&](const std::size_t left, const std::size_t right) {
        for (std::size_t i = left; i < right; ++i) {
            const u64 v = out.index.values[i];
            gcube[i] = sub_mod(sum_cubes_1_to_mod(v), 1ULL);      // sum_{k=2..v} k^3
            gchi[i] = sub_mod(prefix_chi_1_to_mod(v), 1ULL);      // sum_{k=2..v} chi(k)
        }
    });

    std::size_t limit = m;
    for (const int p_int : primes) {
        const u64 p = static_cast<u64>(p_int);
        const u64 p2 = p * p;
        if (p2 > n) {
            break;
        }

        while (limit > 0U && out.index.values[limit - 1U] < p2) {
            --limit;
        }
        if (limit == 0U) {
            break;
        }

        const int ipm1 = out.index.idx_small[static_cast<std::size_t>(p - 1ULL)];
        const u64 gc_pm1 = gcube[static_cast<std::size_t>(ipm1)];
        const u64 gh_pm1 = gchi[static_cast<std::size_t>(ipm1)];

        const u64 p_mod = p % kMod;
        const u64 p3 = mul_mod(mul_mod(p_mod, p_mod), p_mod);
        const int chip = chi4(p);

        auto update_range = [&](const std::size_t left, const std::size_t right) {
            for (std::size_t i = left; i < right; ++i) {
                const u64 v = out.index.values[i];
                const u64 u = v / p;
                const int iu = (u <= out.index.root)
                                   ? out.index.idx_small[static_cast<std::size_t>(u)]
                                   : out.index.idx_large[static_cast<std::size_t>(n / u)];

                const u64 delta_cube = sub_mod(gcube[static_cast<std::size_t>(iu)], gc_pm1);
                gcube[i] = sub_mod(gcube[i], mul_mod(p3, delta_cube));

                if (chip != 0) {
                    const u64 delta_chi = sub_mod(gchi[static_cast<std::size_t>(iu)], gh_pm1);
                    if (chip > 0) {
                        gchi[i] = sub_mod(gchi[i], delta_chi);
                    } else {
                        gchi[i] = add_mod(gchi[i], delta_chi);
                    }
                }
            }
        };

        update_range(0U, limit);
    }

    out.gprime.assign(m, 0ULL);
    parallel_chunked_for(threads, 0U, m, [&](const std::size_t left, const std::size_t right) {
        for (std::size_t i = left; i < right; ++i) {
            out.gprime[i] = sub_mod(gcube[i], gchi[i]);
        }
    });

    return out;
}

class Euler715Solver {
public:
    Euler715Solver(const u64 n, const unsigned threads)
        : n_(n),
          root_(isqrt_u64(n_)),
          primes_(sieve_primes(static_cast<int>(root_))) {
        const PrimeSumData sums = compute_prime_sums_gprime(n_, primes_, threads);
        index_ = sums.index;
        gprime_ = sums.gprime;

        p_count_ = static_cast<int>(primes_.size());
        base_ = static_cast<u64>(p_count_ + 1);
        memo_.reserve(1 << 20);
    }

    u64 solve() {
        return summatory_with_prime_floor(n_, 0);
    }

private:
    u64 n_ = 0ULL;
    u64 root_ = 0ULL;
    std::vector<int> primes_;
    ValuesIndex index_;
    std::vector<u64> gprime_;
    int p_count_ = 0;
    u64 base_ = 1ULL;
    std::unordered_map<u64, u64, U64Hash> memo_;

    int index_of(const u64 x) const {
        if (x <= root_) {
            return index_.idx_small[static_cast<std::size_t>(x)];
        }
        return index_.idx_large[static_cast<std::size_t>(n_ / x)];
    }

    u64 prime_sum_upto(const u64 x) const {
        if (x < 2ULL) {
            return 0ULL;
        }
        return gprime_[static_cast<std::size_t>(index_of(x))];
    }

    u64 prime_sum_range(const u64 lo, const u64 hi) const {
        if (hi <= lo) {
            return 0ULL;
        }
        return sub_mod(prime_sum_upto(hi), prime_sum_upto(lo));
    }

    static u64 g_prime(const u64 p) {
        const u64 p_mod = p % kMod;
        const u64 p3 = mul_mod(mul_mod(p_mod, p_mod), p_mod);
        const int ch = chi4(p);
        if (ch == 1) {
            return sub_mod(p3, 1ULL);
        }
        if (ch == -1) {
            return add_mod(p3, 1ULL);
        }
        return p3;
    }

    u64 summatory_with_prime_floor(const u64 n, const int idx) {
        if (n < 2ULL) {
            return 1ULL;
        }

        if (idx >= p_count_) {
            const u64 lo = (p_count_ == 0) ? 1ULL
                                           : static_cast<u64>(
                                                 primes_[static_cast<std::size_t>(p_count_ - 1)]);
            return add_mod(1ULL, prime_sum_range(lo, n));
        }

        const u64 p0 = static_cast<u64>(primes_[static_cast<std::size_t>(idx)]);
        if (p0 > n) {
            return 1ULL;
        }

        const u64 key = n * base_ + static_cast<u64>(idx);
        const auto it = memo_.find(key);
        if (it != memo_.end()) {
            return it->second;
        }

        const u64 lo = (idx > 0) ? static_cast<u64>(primes_[static_cast<std::size_t>(idx - 1)]) : 1ULL;
        if (p0 > n / p0) {
            const u64 res = add_mod(1ULL, prime_sum_range(lo, n));
            memo_.emplace(key, res);
            return res;
        }

        u64 result = add_mod(1ULL, prime_sum_range(lo, n));

        for (int j = idx; j < p_count_; ++j) {
            const u64 p = static_cast<u64>(primes_[static_cast<std::size_t>(j)]);
            const u64 p2 = p * p;
            if (p2 > n) {
                break;
            }

            const u64 gp1 = g_prime(p);
            const u64 sub = sub_mod(summatory_with_prime_floor(n / p, j + 1), 1ULL);
            result = add_mod(result, mul_mod(gp1, sub));

            const int sign = chi4(p);
            const u64 p_mod = p % kMod;
            const u64 p3 = mul_mod(mul_mod(p_mod, p_mod), p_mod);

            u64 prev = p3;              // p^(3*1)
            u64 cur = mul_mod(prev, p3);  // p^(3*2)
            u64 pe = p2;

            while (pe <= n) {
                u64 gp = 0ULL;
                if (sign == 1) {
                    gp = sub_mod(cur, prev);
                } else if (sign == -1) {
                    gp = add_mod(cur, prev);
                } else {
                    gp = cur;
                }

                result = add_mod(result, mul_mod(gp, summatory_with_prime_floor(n / pe, j + 1)));

                if (pe > n / p) {
                    break;
                }
                pe *= p;
                prev = cur;
                cur = mul_mod(cur, p3);
            }
        }

        memo_.emplace(key, result);
        return result;
    }
};

bool run_checkpoints(const unsigned threads) {
    {
        Euler715Solver solver(kCheckpointN1, threads);
        const u64 got = solver.solve();
        if (got != kCheckpointExpected1) {
            std::cerr << "Validation failed: G(10) expected " << kCheckpointExpected1
                      << ", got " << got << '\n';
            return false;
        }
    }

    {
        Euler715Solver solver(kCheckpointN2, threads);
        const u64 got = solver.solve();
        if (got != kCheckpointExpected2) {
            std::cerr << "Validation failed: G(10^5) expected " << kCheckpointExpected2
                      << ", got " << got << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const unsigned threads = pick_thread_count(options.allow_multithreading, options.requested_threads);

    if (options.run_checkpoints) {
        if (!run_checkpoints(threads)) {
            return 1;
        }
    }

    Euler715Solver solver(options.n, threads);
    const u64 answer = solver.solve();
    std::cout << answer << '\n';
    return 0;
}
