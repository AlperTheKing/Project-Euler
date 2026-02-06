#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultN = 1'000'000'000'000ULL;
constexpr u64 kExpectedAnswerForDefaultN = 132'314'136'838'185ULL;
constexpr u64 kCheckpointN = 1'000'000ULL;
constexpr u64 kCheckpointExpected = 37'429'395ULL;

struct Options {
    u64 n = kDefaultN;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
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

        u64 parsed_u64 = 0ULL;
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            options.n = parsed_u64;
            continue;
        }

        unsigned parsed_unsigned = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_unsigned)) {
            options.requested_threads = parsed_unsigned;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

u64 isqrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) > 0ULL && (r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 icbrt_u64(const u64 x) {
    u64 r = static_cast<u64>(std::cbrt(static_cast<long double>(x)));
    while (static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) *
               static_cast<u128>(r + 1ULL) <=
           static_cast<u128>(x)) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) * static_cast<u128>(r) >
           static_cast<u128>(x)) {
        --r;
    }
    return r;
}

u64 iroot4_u64(const u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(isqrt_u64(x))));
    while (static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) *
               static_cast<u128>(r + 1ULL) * static_cast<u128>(r + 1ULL) <=
           static_cast<u128>(x)) {
        ++r;
    }
    while (static_cast<u128>(r) * static_cast<u128>(r) * static_cast<u128>(r) *
               static_cast<u128>(r) >
           static_cast<u128>(x)) {
        --r;
    }
    return r;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }

    std::string digits;
    while (value > 0U) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(digits.begin(), digits.end());
    return digits;
}

struct PrimeData {
    std::vector<int> primes;
    std::vector<int> pi;

    explicit PrimeData(const int limit) {
        std::vector<int> lp(static_cast<std::size_t>(limit) + 1ULL, 0);
        pi.assign(static_cast<std::size_t>(limit) + 1ULL, 0);

        primes.reserve(static_cast<std::size_t>(limit / 10));
        for (int i = 2; i <= limit; ++i) {
            if (lp[static_cast<std::size_t>(i)] == 0) {
                lp[static_cast<std::size_t>(i)] = i;
                primes.push_back(i);
            }
            for (const int p : primes) {
                const u64 v = static_cast<u64>(p) * static_cast<u64>(i);
                if (v > static_cast<u64>(limit) || p > lp[static_cast<std::size_t>(i)]) {
                    break;
                }
                lp[static_cast<std::size_t>(v)] = p;
            }
            pi[static_cast<std::size_t>(i)] =
                pi[static_cast<std::size_t>(i - 1)] + (lp[static_cast<std::size_t>(i)] == i ? 1 : 0);
        }
    }
};

class LehmerPrimeCounter {
public:
    explicit LehmerPrimeCounter(const PrimeData& data) : data_(data) {}

    u64 pi_leq(const u64 x) {
        if (x < data_.pi.size()) {
            return static_cast<u64>(data_.pi[static_cast<std::size_t>(x)]);
        }

        const auto it = cache_pi_.find(x);
        if (it != cache_pi_.end()) {
            return it->second;
        }

        const u64 a = pi_leq(iroot4_u64(x));
        const u64 b = pi_leq(isqrt_u64(x));
        const u64 c = pi_leq(icbrt_u64(x));

        u64 sum = phi(x, static_cast<int>(a)) + (b + a - 2ULL) * (b - a + 1ULL) / 2ULL;

        for (u64 i = a + 1ULL; i <= b; ++i) {
            const u64 w = x / static_cast<u64>(data_.primes[static_cast<std::size_t>(i - 1ULL)]);
            sum -= pi_leq(w);
            if (i <= c) {
                const u64 bi = pi_leq(isqrt_u64(w));
                for (u64 j = i; j <= bi; ++j) {
                    sum -= pi_leq(w / static_cast<u64>(data_.primes[static_cast<std::size_t>(j - 1ULL)])) -
                           (j - 1ULL);
                }
            }
        }

        cache_pi_.emplace(x, sum);
        return sum;
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
            return x - x / 2ULL - x / 3ULL - x / 5ULL + x / 6ULL + x / 10ULL + x / 15ULL -
                   x / 30ULL;
        }

        if (x < data_.pi.size() &&
            static_cast<u64>(data_.primes[static_cast<std::size_t>(s - 1)]) *
                    static_cast<u64>(data_.primes[static_cast<std::size_t>(s - 1)]) >
                x) {
            return static_cast<u64>(data_.pi[static_cast<std::size_t>(x)]) - static_cast<u64>(s) +
                   1ULL;
        }

        // Tiny-s cache keeps Lehmer fast while avoiding a huge memo footprint.
        if (s <= 6 && x <= 1'000'000'000ULL) {
            const u64 key = (x << 6U) ^ static_cast<u64>(s);
            const auto it = cache_phi_.find(key);
            if (it != cache_phi_.end()) {
                return it->second;
            }
            const u64 result =
                phi(x, s - 1) - phi(x / static_cast<u64>(data_.primes[static_cast<std::size_t>(s - 1)]),
                                    s - 1);
            cache_phi_.emplace(key, result);
            return result;
        }

        return phi(x, s - 1) -
               phi(x / static_cast<u64>(data_.primes[static_cast<std::size_t>(s - 1)]), s - 1);
    }

    const PrimeData& data_;
    std::unordered_map<u64, u64> cache_phi_;
    std::unordered_map<u64, u64> cache_pi_;
};

class SummatoryD2Solver {
public:
    explicit SummatoryD2Solver(const PrimeData& data)
        : data_(data), prime_counter_(data) {}

    std::size_t root_branch_count(const u64 n) const {
        const u64 root = isqrt_u64(n);
        return static_cast<std::size_t>(std::upper_bound(data_.primes.begin(),
                                                         data_.primes.end(),
                                                         static_cast<int>(root)) -
                                        data_.primes.begin());
    }

    u128 solve_single(const u64 n) {
        return restricted_sum(n, 0) + 1U;
    }

    u128 solve_parallel_root(const u64 n, const unsigned threads) {
        if (threads <= 1U) {
            return solve_single(n);
        }

        SummatoryD2Solver base_solver(data_);
        u128 total = 1U + static_cast<u128>(3ULL) * static_cast<u128>(base_solver.prime_counter_.pi_leq(n));

        const int root_count = static_cast<int>(root_branch_count(n));
        if (root_count == 0) {
            return total;
        }

        const unsigned use_threads = std::min<unsigned>(threads, static_cast<unsigned>(root_count));
        if (use_threads == 1U) {
            return solve_single(n);
        }

        std::vector<u128> partial(static_cast<std::size_t>(use_threads), 0U);
        std::vector<std::thread> pool;
        pool.reserve(static_cast<std::size_t>(use_threads));

        for (unsigned t = 0U; t < use_threads; ++t) {
            pool.emplace_back([&, t]() {
                SummatoryD2Solver worker(data_);
                u128 local = 0U;
                for (int i = static_cast<int>(t); i < root_count; i += static_cast<int>(use_threads)) {
                    local += worker.root_branch_contribution(n, i);
                }
                partial[static_cast<std::size_t>(t)] = local;
            });
        }

        for (std::thread& th : pool) {
            th.join();
        }

        for (const u128 part : partial) {
            total += part;
        }
        return total;
    }

private:
    u128 root_branch_contribution(const u64 n, const int prime_index) {
        const u64 p = static_cast<u64>(data_.primes[static_cast<std::size_t>(prime_index)]);
        u64 power = p;
        int exponent = 1;
        u128 sum = 0U;

        while (power <= n) {
            u128 tail = restricted_sum(n / power, prime_index + 1);
            if (exponent >= 2) {
                tail += 1U;
            }
            sum += static_cast<u128>(2 * exponent + 1) * tail;

            if (power > n / p) {
                break;
            }
            power *= p;
            ++exponent;
        }

        return sum;
    }

    u128 restricted_sum(const u64 n, const int min_prime_index) {
        if (n < 2ULL) {
            return 0U;
        }
        if (min_prime_index < static_cast<int>(data_.primes.size()) &&
            static_cast<u64>(data_.primes[static_cast<std::size_t>(min_prime_index)]) > n) {
            return 0U;
        }

        u128 answer = 0U;
        const u64 prime_count = prime_counter_.pi_leq(n);
        if (prime_count > static_cast<u64>(min_prime_index)) {
            answer += static_cast<u128>(3ULL) *
                      static_cast<u128>(prime_count - static_cast<u64>(min_prime_index));
        }

        for (int i = min_prime_index; i < static_cast<int>(data_.primes.size()); ++i) {
            const u64 p = static_cast<u64>(data_.primes[static_cast<std::size_t>(i)]);
            if (p * p > n) {
                break;
            }

            u64 power = p;
            int exponent = 1;

            while (power <= n) {
                u128 tail = restricted_sum(n / power, i + 1);
                if (exponent >= 2) {
                    tail += 1U;
                }
                answer += static_cast<u128>(2 * exponent + 1) * tail;

                if (power > n / p) {
                    break;
                }
                power *= p;
                ++exponent;
            }
        }

        return answer;
    }

    const PrimeData& data_;
    LehmerPrimeCounter prime_counter_;
};

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    constexpr std::size_t kMinUnitsForParallel = 64ULL;
    constexpr unsigned kDefaultThreadCap = 8U;

    if (!allow_multithreading || workload_units < 2ULL || workload_units < kMinUnitsForParallel) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
        threads = std::min(threads, kDefaultThreadCap);
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

u64 gcd_u64(u64 a, u64 b) {
    while (b != 0ULL) {
        const u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u64 brute_f_from_lcm(const u64 n) {
    u64 count = 0ULL;
    for (u64 x = 1ULL; x <= n; ++x) {
        for (u64 y = x; y <= n; ++y) {
            const u64 g = gcd_u64(x, y);
            const u64 l = (x / g) * y;
            if (l == n) {
                ++count;
            }
        }
    }
    return count;
}

u64 formula_f_from_factorization(u64 n) {
    u64 divisor_count_n2 = 1ULL;
    for (u64 p = 2ULL; p * p <= n; ++p) {
        if (n % p != 0ULL) {
            continue;
        }

        int exponent = 0;
        while (n % p == 0ULL) {
            n /= p;
            ++exponent;
        }
        divisor_count_n2 *= static_cast<u64>(2 * exponent + 1);
    }
    if (n > 1ULL) {
        divisor_count_n2 *= 3ULL;
    }
    return (divisor_count_n2 + 1ULL) / 2ULL;
}

u64 summatory_g_via_spf(const int n) {
    std::vector<int> spf(static_cast<std::size_t>(n) + 1ULL, 0);
    for (int i = 2; i <= n; ++i) {
        if (spf[static_cast<std::size_t>(i)] != 0) {
            continue;
        }
        for (u64 j = static_cast<u64>(i); j <= static_cast<u64>(n); j += static_cast<u64>(i)) {
            if (spf[static_cast<std::size_t>(j)] == 0) {
                spf[static_cast<std::size_t>(j)] = i;
            }
        }
    }

    u64 total = 0ULL;
    for (int x = 1; x <= n; ++x) {
        int value = x;
        u64 divisor_count_n2 = 1ULL;
        while (value > 1) {
            const int p = spf[static_cast<std::size_t>(value)];
            int exponent = 0;
            while (value % p == 0) {
                value /= p;
                ++exponent;
            }
            divisor_count_n2 *= static_cast<u64>(2 * exponent + 1);
        }
        total += (divisor_count_n2 + 1ULL) / 2ULL;
    }
    return total;
}

bool run_checkpoints(const PrimeData& prime_data) {
    for (u64 n = 1ULL; n <= 200ULL; ++n) {
        const u64 brute = brute_f_from_lcm(n);
        const u64 formula = formula_f_from_factorization(n);
        if (brute != formula) {
            std::cerr << "Identity check failed at n=" << n << ": brute f(n)=" << brute
                      << ", formula f(n)=" << formula << '\n';
            return false;
        }
    }

    struct Checkpoint {
        int n;
        u64 expected;
    };
    const std::vector<Checkpoint> small = {
        {10, 29ULL},
        {100, 647ULL},
        {1'000, 11'751ULL},
        {10'000, 186'991ULL},
    };

    for (const Checkpoint& cp : small) {
        const u64 got = summatory_g_via_spf(cp.n);
        if (got != cp.expected) {
            std::cerr << "Small checkpoint failed for g(" << cp.n << "): expected " << cp.expected
                      << ", got " << got << '\n';
            return false;
        }
    }

    SummatoryD2Solver single_solver(prime_data);
    const u128 single_s = single_solver.solve_single(kCheckpointN);
    const u128 single_g = (single_s + static_cast<u128>(kCheckpointN)) / 2U;
    if (static_cast<u64>(single_g) != kCheckpointExpected) {
        std::cerr << "Main checkpoint failed for g(" << kCheckpointN << "): expected "
                  << kCheckpointExpected << ", got " << to_string_u128(single_g) << '\n';
        return false;
    }

    SummatoryD2Solver parallel_solver(prime_data);
    const u128 parallel_s = parallel_solver.solve_parallel_root(kCheckpointN, 2U);
    if (parallel_s != single_s) {
        std::cerr << "Parallel consistency check failed at n=" << kCheckpointN << '\n';
        return false;
    }

    return true;
}

int choose_prime_limit(const u64 n) {
    const u64 need = std::max<u64>(5'000'000ULL, isqrt_u64(n) + 100ULL);
    if (need > static_cast<u64>(std::numeric_limits<int>::max())) {
        return -1;
    }
    return static_cast<int>(need);
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const int prime_limit = choose_prime_limit(options.n);
    if (prime_limit < 0) {
        std::cerr << "Requested n is too large for this implementation's prime table.\n";
        return 1;
    }

    const PrimeData prime_data(prime_limit);

    if (options.run_checkpoints && !run_checkpoints(prime_data)) {
        return 1;
    }

    SummatoryD2Solver solver(prime_data);
    const std::size_t workload = solver.root_branch_count(options.n);
    const unsigned threads =
        choose_thread_count(options.allow_multithreading, options.requested_threads, workload);

    const u128 summatory_d_n_square =
        (threads == 1U) ? solver.solve_single(options.n)
                        : solver.solve_parallel_root(options.n, threads);
    const u128 answer = (summatory_d_n_square + static_cast<u128>(options.n)) / 2U;

    if (options.n == kDefaultN && static_cast<u64>(answer) != kExpectedAnswerForDefaultN) {
        std::cerr << "Internal validation failed for n=" << kDefaultN << ": expected "
                  << kExpectedAnswerForDefaultN << ", got " << to_string_u128(answer) << '\n';
        return 1;
    }

    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
