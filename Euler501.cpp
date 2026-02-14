#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct Options {
    u64 n = 1'000'000'000'000ULL;
    bool run_checkpoints = true;
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
    for (char ch : tail) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        parsed = parsed * 10ULL + static_cast<u64>(ch - '0');
    }
    value = parsed;
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
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
    return options.n >= 1ULL;
}

u64 isqrt_u64(const u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / (x + 1ULL)) {
        ++x;
    }
    while (x > 0ULL && x > n / x) {
        --x;
    }
    return x;
}

u64 icbrt_u64(const u64 n) {
    u64 x = static_cast<u64>(std::cbrt(static_cast<long double>(n)));
    while ((x + 1ULL) <= n / ((x + 1ULL) * (x + 1ULL))) {
        ++x;
    }
    while (x > 0ULL && x > n / (x * x)) {
        --x;
    }
    return x;
}

u64 iroot4_u64(const u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(std::sqrt(static_cast<long double>(n)))));
    auto pow4 = [](u64 y) -> __uint128_t {
        return static_cast<__uint128_t>(y) * y * y * y;
    };
    while (pow4(x + 1ULL) <= n) {
        ++x;
    }
    while (x > 0ULL && pow4(x) > n) {
        --x;
    }
    return x;
}

class PrimeCounting {
public:
    PrimeCounting() {
        constexpr int MAX = 5'000'000;
        sieve_limit_ = MAX;
        std::vector<bool> is_comp(static_cast<std::size_t>(MAX + 1), false);
        pi_small_.assign(static_cast<std::size_t>(MAX + 1), 0);
        for (int i = 2; i <= MAX; ++i) {
            if (!is_comp[static_cast<std::size_t>(i)]) {
                primes_.push_back(i);
                if (i <= MAX / i) {
                    for (int j = i * i; j <= MAX; j += i) {
                        is_comp[static_cast<std::size_t>(j)] = true;
                    }
                }
            }
            pi_small_[static_cast<std::size_t>(i)] = pi_small_[static_cast<std::size_t>(i - 1)] +
                                                     (!is_comp[static_cast<std::size_t>(i)] ? 1 : 0);
        }
        for (int i = 0; i < 7; ++i) {
            primes_for_phi_.push_back(primes_[static_cast<std::size_t>(i)]);
        }
    }

    u64 pi(const u64 n) {
        if (n <= static_cast<u64>(sieve_limit_)) {
            return pi_small_[static_cast<std::size_t>(n)];
        }
        auto it = pi_cache_.find(n);
        if (it != pi_cache_.end()) {
            return it->second;
        }

        const u64 a = pi(iroot4_u64(n));
        const u64 b = pi(isqrt_u64(n));
        const u64 c = pi(icbrt_u64(n));

        i64 sum = static_cast<i64>(phi(n, static_cast<int>(a))) +
                  static_cast<i64>((b + a - 2ULL) * (b - a + 1ULL) / 2ULL);

        for (u64 i = a + 1ULL; i <= b; ++i) {
            const u64 w = n / static_cast<u64>(primes_[static_cast<std::size_t>(i - 1ULL)]);
            sum -= static_cast<i64>(pi(w));
            if (i <= c) {
                const u64 lim = pi(isqrt_u64(w));
                for (u64 j = i; j <= lim; ++j) {
                    const u64 pj = static_cast<u64>(primes_[static_cast<std::size_t>(j - 1ULL)]);
                    sum -= static_cast<i64>(pi(w / pj) - (j - 1ULL));
                }
            }
        }

        const u64 out = static_cast<u64>(sum);
        pi_cache_[n] = out;
        return out;
    }

    const std::vector<int>& primes() const {
        return primes_;
    }

private:
    u64 phi(const u64 x, const int s) {
        if (s == 0) {
            return x;
        }
        if (s == 1) {
            return x - x / 2ULL;
        }
        if (s == 2) {
            return x - x / 2ULL - x / 3ULL + x / 6ULL;
        }
        if (s == 3) {
            return x - x / 2ULL - x / 3ULL - x / 5ULL + x / 6ULL + x / 10ULL + x / 15ULL - x / 30ULL;
        }
        if (x <= static_cast<u64>(sieve_limit_) && static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]) >= x) {
            return 1ULL;
        }
        const u64 key = (x << 16U) ^ static_cast<u64>(s);
        auto it = phi_cache_.find(key);
        if (it != phi_cache_.end()) {
            return it->second;
        }
        const u64 out = phi(x, s - 1) - phi(x / static_cast<u64>(primes_[static_cast<std::size_t>(s - 1)]), s - 1);
        phi_cache_[key] = out;
        return out;
    }

    int sieve_limit_ = 0;
    std::vector<int> primes_;
    std::vector<u64> pi_small_;
    std::vector<int> primes_for_phi_;
    std::unordered_map<u64, u64> pi_cache_;
    std::unordered_map<u64, u64> phi_cache_;
};

u64 solve(const u64 n) {
    PrimeCounting pc;
    const std::vector<int>& primes = pc.primes();

    u64 total = 0ULL;

    // Form p^7
    total += pc.pi(static_cast<u64>(std::pow(static_cast<long double>(n), 1.0L / 7.0L)));
    while (true) {
        const u64 p = pc.pi(total);
        (void)p;
        break;
    }
    u64 p7_count = 0ULL;
    for (int p : primes) {
        __uint128_t p7 = 1;
        for (int i = 0; i < 7; ++i) {
            p7 *= static_cast<u64>(p);
            if (p7 > n) {
                break;
            }
        }
        if (p7 <= n) {
            ++p7_count;
        } else {
            break;
        }
    }
    total = p7_count;

    // Form p^3 q
    for (int p : primes) {
        const u64 p64 = static_cast<u64>(p);
        const u64 p3 = p64 * p64 * p64;
        if (p3 > n) {
            break;
        }
        const u64 qmax = n / p3;
        if (qmax < 2ULL) {
            continue;
        }
        u64 cnt = pc.pi(qmax);
        if (p64 <= qmax) {
            --cnt;
        }
        total += cnt;
    }

    // Form p q r, distinct primes, p<q<r
    const std::size_t psz = primes.size();
    for (std::size_t i = 0; i < psz; ++i) {
        const u64 p = static_cast<u64>(primes[i]);
        if (p * p * p > n) {
            break;
        }
        for (std::size_t j = i + 1; j < psz; ++j) {
            const u64 q = static_cast<u64>(primes[j]);
            if (p > n / (q * q)) {
                break;
            }
            const u64 rmax = n / (p * q);
            if (rmax <= q) {
                continue;
            }
            const u64 upto = pc.pi(rmax);
            const u64 before_or_equal_q = static_cast<u64>(j + 1);
            if (upto > before_or_equal_q) {
                total += (upto - before_or_equal_q);
            }
        }
    }

    return total;
}

u64 brute_count(const u64 n) {
    u64 count = 0ULL;
    for (u64 x = 1ULL; x <= n; ++x) {
        u64 y = x;
        int d = 1;
        for (u64 p = 2ULL; p * p <= y; ++p) {
            if (y % p != 0ULL) {
                continue;
            }
            int e = 0;
            while (y % p == 0ULL) {
                y /= p;
                ++e;
            }
            d *= (e + 1);
        }
        if (y > 1ULL) {
            d *= 2;
        }
        if (d == 8) {
            ++count;
        }
    }
    return count;
}

bool run_checkpoints() {
    if (solve(100ULL) != 10ULL) {
        std::cerr << "Checkpoint failed: f(100)=10" << '\n';
        return false;
    }
    if (solve(1000ULL) != 180ULL) {
        std::cerr << "Checkpoint failed: f(1000)=180" << '\n';
        return false;
    }
    if (solve(1'000'000ULL) != 224'427ULL) {
        std::cerr << "Checkpoint failed: f(1e6)=224427" << '\n';
        return false;
    }
    if (solve(20'000ULL) != brute_count(20'000ULL)) {
        std::cerr << "Checkpoint failed: brute-force cross-check for n=20000" << '\n';
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
        return 2;
    }
    std::cout << solve(options.n) << '\n';
    return 0;
}
