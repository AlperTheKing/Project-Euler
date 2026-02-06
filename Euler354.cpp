#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTargetL = 500'000'000'000ULL;
constexpr u64 kTargetB = 450ULL;
constexpr u64 kTargetDiv6 = kTargetB / 6ULL;  // 75
constexpr u128 kMinCore = static_cast<u128>(2'401ULL) * 28'561ULL * 361ULL;       // 7^4 * 13^4 * 19^2
constexpr u128 kCoreForExp2Bound = static_cast<u128>(2'401ULL) * 28'561ULL;        // 7^4 * 13^4

struct Options {
    u64 limit_l = kTargetL;
    bool run_checks = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) return false;

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) return false;

    u64 parsed = 0ULL;
    for (const char c : tail) {
        if (c < '0' || c > '9') return false;
        const u64 digit = static_cast<u64>(c - '0');
        if (parsed > (std::numeric_limits<u64>::max() - digit) / 10ULL) return false;
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) return false;
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) return false;
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checks") {
            options.run_checks = false;
            continue;
        }
        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }

        unsigned parsed_threads = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_threads)) {
            options.requested_threads = parsed_threads;
            continue;
        }

        u64 parsed_l = 0ULL;
        if (parse_u64_after_prefix(arg, "--limit-l=", parsed_l)) {
            options.limit_l = parsed_l;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.limit_l == 0ULL) {
        std::cerr << "L must be positive.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    if (!allow_multithreading || workload_units < 2ULL) return 1U;

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) threads = 1U;
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

u128 from_u64(const u64 x) {
    return static_cast<u128>(x);
}

std::string to_string_u128(u128 value) {
    if (value == 0) return "0";

    std::string out;
    while (value > 0) {
        const u32 digit = static_cast<u32>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

u128 pow_u128_limited(const u64 base, const int exponent, const u128 limit) {
    u128 out = 1;
    const u128 b = from_u64(base);
    for (int i = 0; i < exponent; ++i) {
        if (out > limit / b) return limit + 1;
        out *= b;
    }
    return out;
}

u64 isqrt_u128(const u128 value) {
    if (value == 0) return 0ULL;

    long double approx = std::sqrt(static_cast<long double>(value));
    u64 root = static_cast<u64>(approx);
    while (from_u64(root + 1ULL) * from_u64(root + 1ULL) <= value) ++root;
    while (from_u64(root) * from_u64(root) > value) --root;
    return root;
}

u128 square_u64(const u64 x) {
    return from_u64(x) * from_u64(x);
}

std::vector<int> sieve_primes_upto(const int limit) {
    std::vector<int> primes;
    if (limit < 2) return primes;

    std::vector<std::uint8_t> is_composite(static_cast<std::size_t>(limit + 1), 0U);
    primes.reserve(static_cast<std::size_t>(limit / 10));

    for (int i = 2; i <= limit; ++i) {
        if (is_composite[static_cast<std::size_t>(i)] == 0U) {
            primes.push_back(i);
            if (from_u64(static_cast<u64>(i)) * from_u64(static_cast<u64>(i)) <= from_u64(static_cast<u64>(limit))) {
                for (u64 j = from_u64(static_cast<u64>(i)) * from_u64(static_cast<u64>(i));
                     j <= static_cast<u64>(limit);
                     j += static_cast<u64>(i)) {
                    is_composite[static_cast<std::size_t>(j)] = 1U;
                }
            }
        }
    }

    return primes;
}

std::vector<int> sieve_primes_mod_1_upto(const int limit) {
    const std::vector<int> all_primes = sieve_primes_upto(limit);
    std::vector<int> mod1;
    mod1.reserve(all_primes.size() / 2U + 8U);
    for (const int p : all_primes) {
        if (p % 3 == 1) mod1.push_back(p);
    }
    return mod1;
}

u64 count_representations_formula(u64 n, const std::vector<int>& factor_primes) {
    u64 product = 1ULL;

    for (const int prime : factor_primes) {
        const u64 p = static_cast<u64>(prime);
        if (p * p > n) break;
        if (n % p != 0ULL) continue;

        int exponent = 0;
        while (n % p == 0ULL) {
            n /= p;
            ++exponent;
        }

        if (p % 3ULL == 1ULL) {
            product *= static_cast<u64>(exponent + 1);
        } else if (p % 3ULL == 2ULL) {
            if ((exponent & 1) != 0) return 0ULL;
        }
    }

    if (n > 1ULL) {
        if (n % 3ULL == 1ULL) {
            product *= 2ULL;
        } else if (n % 3ULL == 2ULL) {
            return 0ULL;
        }
    }

    return 6ULL * product;
}

u64 count_representations_bruteforce(const int target_n) {
    const int radius = static_cast<int>(std::sqrt(static_cast<double>(target_n))) + 2;
    u64 count = 0ULL;
    for (int a = -radius; a <= radius; ++a) {
        for (int b = -radius; b <= radius; ++b) {
            const int value = a * a + a * b + b * b;
            if (value == target_n) ++count;
        }
    }
    return count;
}

bool is_multiplier_bruteforce(u64 value, const std::vector<int>& factor_primes) {
    u64 n = value;
    for (const int prime : factor_primes) {
        const u64 p = static_cast<u64>(prime);
        if (p * p > n) break;
        if (n % p != 0ULL) continue;

        int exponent = 0;
        while (n % p == 0ULL) {
            n /= p;
            ++exponent;
        }

        if (p == 3ULL) continue;
        if (p % 3ULL == 1ULL) return false;
        if ((exponent & 1) != 0) return false;
    }

    if (n > 1ULL) {
        if (n == 3ULL) return true;
        if (n % 3ULL == 1ULL) return false;
        return false;  // n % 3 == 2 with odd exponent 1.
    }

    return true;
}

struct NeutralCounter {
    std::vector<u32> prefix;

    u64 count_multipliers(const u128 x) const {
        u64 total = 0ULL;
        u128 power3 = 1;

        while (power3 <= x) {
            const u64 y = isqrt_u128(x / power3);
            total += prefix[static_cast<std::size_t>(y)];
            if (power3 > x / 3U) break;
            power3 *= 3U;
        }

        return total;
    }
};

NeutralCounter build_neutral_counter(const u64 max_s) {
    NeutralCounter out;
    out.prefix.assign(static_cast<std::size_t>(max_s + 1U), 0U);
    if (max_s == 0ULL) return out;

    std::vector<int> spf(static_cast<std::size_t>(max_s + 1U), 0);
    for (u64 i = 2ULL; i <= max_s; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = static_cast<int>(i);
            if (i <= max_s / i) {
                for (u64 j = i * i; j <= max_s; j += i) {
                    if (spf[static_cast<std::size_t>(j)] == 0) {
                        spf[static_cast<std::size_t>(j)] = static_cast<int>(i);
                    }
                }
            }
        }
    }

    std::vector<u32> good(static_cast<std::size_t>(max_s + 1U), 0U);
    good[1] = 1U;

    for (u64 n = 2ULL; n <= max_s; ++n) {
        const int p = spf[static_cast<std::size_t>(n)];
        const u64 m = n / static_cast<u64>(p);
        good[static_cast<std::size_t>(n)] =
            (p % 3 == 2 && good[static_cast<std::size_t>(m)] == 1U) ? 1U : 0U;
    }

    for (u64 n = 1ULL; n <= max_s; ++n) {
        out.prefix[static_cast<std::size_t>(n)] =
            out.prefix[static_cast<std::size_t>(n - 1U)] + good[static_cast<std::size_t>(n)];
    }

    return out;
}

template <typename Fn>
u128 parallel_sum_indices(const std::size_t i_end,
                          const bool allow_multithreading,
                          const unsigned requested_threads,
                          Fn&& fn) {
    if (i_end == 0U) return 0;

    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, i_end);
    std::vector<u128> partial(threads, 0);
    std::atomic<std::size_t> next_i{0U};

    const auto worker = [&](const unsigned tid) {
        u128 local = 0;
        while (true) {
            const std::size_t i = next_i.fetch_add(1U, std::memory_order_relaxed);
            if (i >= i_end) break;
            local += fn(i);
        }
        partial[tid] = local;
    };

    if (threads == 1U) {
        worker(0U);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(threads);
        for (unsigned t = 0U; t < threads; ++t) {
            pool.emplace_back(worker, t);
        }
        for (std::thread& th : pool) th.join();
    }

    u128 total = 0;
    for (const u128 value : partial) total += value;
    return total;
}

u128 solve_for_n_limit(const u128 n_limit,
                       const bool allow_multithreading,
                       const unsigned requested_threads) {
    if (n_limit < kMinCore) return 0;
    const int prime_limit = static_cast<int>(isqrt_u128(n_limit / kCoreForExp2Bound));

    const std::vector<int> primes_mod1 = sieve_primes_mod_1_upto(prime_limit);
    const std::size_t ps = primes_mod1.size();

    std::vector<u128> pow2(ps, n_limit + 1);
    std::vector<u128> pow4(ps, n_limit + 1);
    for (std::size_t i = 0; i < ps; ++i) {
        const u64 p = static_cast<u64>(primes_mod1[i]);
        pow2[i] = pow_u128_limited(p, 2, n_limit);
        pow4[i] = pow_u128_limited(p, 4, n_limit);
    }

    std::vector<std::pair<std::size_t, u128>> pow14_list;
    std::vector<std::pair<std::size_t, u128>> pow24_list;
    std::vector<std::pair<std::size_t, u128>> pow74_list;

    for (std::size_t i = 0; i < ps; ++i) {
        const u64 p = static_cast<u64>(primes_mod1[i]);
        const u128 p14 = pow_u128_limited(p, 14, n_limit);
        if (p14 <= n_limit) {
            pow14_list.push_back({i, p14});
        } else {
            break;
        }
    }
    for (std::size_t i = 0; i < ps; ++i) {
        const u64 p = static_cast<u64>(primes_mod1[i]);
        const u128 p24 = pow_u128_limited(p, 24, n_limit);
        if (p24 <= n_limit) {
            pow24_list.push_back({i, p24});
        } else {
            break;
        }
    }
    for (std::size_t i = 0; i < ps; ++i) {
        const u64 p = static_cast<u64>(primes_mod1[i]);
        const u128 p74 = pow_u128_limited(p, 74, n_limit);
        if (p74 <= n_limit) {
            pow74_list.push_back({i, p74});
        } else {
            break;
        }
    }

    const u64 max_s = isqrt_u128(n_limit / kMinCore);
    const NeutralCounter neutral = build_neutral_counter(max_s);

    const std::size_t i_end4 = static_cast<std::size_t>(
        std::upper_bound(pow4.begin(), pow4.end(), n_limit) - pow4.begin());

    u128 answer = 0;

    // Pattern [74].
    for (const auto& entry : pow74_list) {
        answer += neutral.count_multipliers(n_limit / entry.second);
    }

    // Pattern [24, 2].
    for (const auto& entry : pow24_list) {
        const std::size_t i = entry.first;
        const u128 b0 = n_limit / entry.second;
        const std::size_t j_begin = i + 1U;
        if (j_begin >= ps) continue;
        const std::size_t j_end = static_cast<std::size_t>(
            std::upper_bound(pow2.begin() + static_cast<std::ptrdiff_t>(j_begin), pow2.end(), b0) -
            pow2.begin());
        for (std::size_t j = j_begin; j < j_end; ++j) {
            answer += neutral.count_multipliers(b0 / pow2[j]);
        }
    }

    // Pattern [2, 24].
    for (const auto& entry : pow24_list) {
        const std::size_t j = entry.first;
        if (j == 0U) continue;
        const u128 b0 = n_limit / entry.second;
        const std::size_t i_end = static_cast<std::size_t>(
            std::upper_bound(pow2.begin(), pow2.begin() + static_cast<std::ptrdiff_t>(j), b0) -
            pow2.begin());
        for (std::size_t i = 0; i < i_end; ++i) {
            answer += neutral.count_multipliers(b0 / pow2[i]);
        }
    }

    // Pattern [14, 4].
    for (const auto& entry : pow14_list) {
        const std::size_t i = entry.first;
        const u128 b0 = n_limit / entry.second;
        const std::size_t j_begin = i + 1U;
        if (j_begin >= ps) continue;
        const std::size_t j_end = static_cast<std::size_t>(
            std::upper_bound(pow4.begin() + static_cast<std::ptrdiff_t>(j_begin), pow4.end(), b0) -
            pow4.begin());
        for (std::size_t j = j_begin; j < j_end; ++j) {
            answer += neutral.count_multipliers(b0 / pow4[j]);
        }
    }

    // Pattern [4, 14].
    for (const auto& entry : pow14_list) {
        const std::size_t j = entry.first;
        if (j == 0U) continue;
        const u128 b0 = n_limit / entry.second;
        const std::size_t i_end = static_cast<std::size_t>(
            std::upper_bound(pow4.begin(), pow4.begin() + static_cast<std::ptrdiff_t>(j), b0) -
            pow4.begin());
        for (std::size_t i = 0; i < i_end; ++i) {
            answer += neutral.count_multipliers(b0 / pow4[i]);
        }
    }

    // Pattern [4, 4, 2].
    answer += parallel_sum_indices(
        i_end4, allow_multithreading, requested_threads, [&](const std::size_t i) -> u128 {
            const u128 b0 = n_limit / pow4[i];
            const std::size_t j_begin = i + 1U;
            if (j_begin >= ps) return 0;
            const std::size_t j_end = static_cast<std::size_t>(
                std::upper_bound(pow4.begin() + static_cast<std::ptrdiff_t>(j_begin), pow4.end(), b0) -
                pow4.begin());

            u128 local = 0;
            for (std::size_t j = j_begin; j < j_end; ++j) {
                const u128 b1 = b0 / pow4[j];
                const std::size_t k_begin = j + 1U;
                if (k_begin >= ps) continue;
                const std::size_t k_end = static_cast<std::size_t>(
                    std::upper_bound(pow2.begin() + static_cast<std::ptrdiff_t>(k_begin),
                                     pow2.end(),
                                     b1) -
                    pow2.begin());
                for (std::size_t k = k_begin; k < k_end; ++k) {
                    local += neutral.count_multipliers(b1 / pow2[k]);
                }
            }
            return local;
        });

    // Pattern [4, 2, 4].
    answer += parallel_sum_indices(
        i_end4, allow_multithreading, requested_threads, [&](const std::size_t i) -> u128 {
            const u128 b0 = n_limit / pow4[i];
            const std::size_t j_begin = i + 1U;
            if (j_begin >= ps) return 0;
            const std::size_t j_end = static_cast<std::size_t>(
                std::upper_bound(pow2.begin() + static_cast<std::ptrdiff_t>(j_begin), pow2.end(), b0) -
                pow2.begin());

            u128 local = 0;
            for (std::size_t j = j_begin; j < j_end; ++j) {
                const u128 b1 = b0 / pow2[j];
                const std::size_t k_begin = j + 1U;
                if (k_begin >= ps) continue;
                const std::size_t k_end = static_cast<std::size_t>(
                    std::upper_bound(pow4.begin() + static_cast<std::ptrdiff_t>(k_begin),
                                     pow4.end(),
                                     b1) -
                    pow4.begin());
                for (std::size_t k = k_begin; k < k_end; ++k) {
                    local += neutral.count_multipliers(b1 / pow4[k]);
                }
            }
            return local;
        });

    // Pattern [2, 4, 4].
    const std::size_t i_end_for_244 = (i_end4 > 0U) ? (i_end4 - 1U) : 0U;
    answer += parallel_sum_indices(
        i_end_for_244, allow_multithreading, requested_threads, [&](const std::size_t i) -> u128 {
            const u128 b0 = n_limit / pow2[i];
            const std::size_t j_begin = i + 1U;
            if (j_begin >= ps) return 0;
            const std::size_t j_end = static_cast<std::size_t>(
                std::upper_bound(pow4.begin() + static_cast<std::ptrdiff_t>(j_begin), pow4.end(), b0) -
                pow4.begin());

            u128 local = 0;
            for (std::size_t j = j_begin; j < j_end; ++j) {
                const u128 b1 = b0 / pow4[j];
                const std::size_t k_begin = j + 1U;
                if (k_begin >= ps) continue;
                const std::size_t k_end = static_cast<std::size_t>(
                    std::upper_bound(pow4.begin() + static_cast<std::ptrdiff_t>(k_begin),
                                     pow4.end(),
                                     b1) -
                    pow4.begin());
                for (std::size_t k = k_begin; k < k_end; ++k) {
                    local += neutral.count_multipliers(b1 / pow4[k]);
                }
            }
            return local;
        });

    return answer;
}

u128 n_limit_from_l_limit(const u64 l_limit) {
    return square_u64(l_limit) / 3U;
}

bool run_validations(const Options& options) {
    const std::vector<int> small_primes = sieve_primes_upto(5'000'000);

    // Checkpoint 1: formula for r(n) against brute lattice counting.
    for (int n = 1; n <= 120; ++n) {
        const u64 brute = count_representations_bruteforce(n);
        const u64 fast = count_representations_formula(static_cast<u64>(n), small_primes);
        if (brute != fast) {
            std::cerr << "Validation failed: representation mismatch at n=" << n
                      << " (brute=" << brute << ", fast=" << fast << ")\n";
            return false;
        }
    }

    // Checkpoint 2: statement examples.
    if (count_representations_formula(1ULL, small_primes) != 6ULL) {
        std::cerr << "Validation failed: B(sqrt(3)) != 6\n";
        return false;
    }
    if (count_representations_formula(7ULL, small_primes) != 12ULL) {
        std::cerr << "Validation failed: B(sqrt(21)) != 12\n";
        return false;
    }

    const u64 l_example = 111'111'111ULL;
    const u64 n_example = static_cast<u64>(n_limit_from_l_limit(l_example));
    if (count_representations_formula(n_example, small_primes) != 54ULL) {
        std::cerr << "Validation failed: B(111111111) != 54\n";
        return false;
    }

    // Checkpoint 3: neutral multiplier counter against brute force.
    {
        const u64 max_s = 2'000ULL;
        const NeutralCounter neutral = build_neutral_counter(max_s);
        const std::vector<u64> tests = {1ULL, 3ULL, 10ULL, 100ULL, 1'000ULL, 10'000ULL};

        for (const u64 x : tests) {
            u64 brute = 0ULL;
            for (u64 t = 1ULL; t <= x; ++t) {
                if (is_multiplier_bruteforce(t, small_primes)) ++brute;
            }

            const u64 fast = neutral.count_multipliers(from_u64(x));
            if (brute != fast) {
                std::cerr << "Validation failed: neutral multiplier mismatch at x=" << x
                          << " (brute=" << brute << ", fast=" << fast << ")\n";
                return false;
            }
        }
    }

    // Checkpoint 4: single-thread vs multi-thread consistency on a reduced limit.
    if (options.allow_multithreading) {
        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0U) hw = 1U;
        if (hw > 1U) {
            const u64 l_test = 50'000'000ULL;
            const u128 n_test = n_limit_from_l_limit(l_test);
            const u128 single = solve_for_n_limit(n_test, false, 1U);
            const u128 multi = solve_for_n_limit(n_test, true, options.requested_threads);
            if (single != multi) {
                std::cerr << "Validation failed: thread consistency mismatch at L=" << l_test
                          << " (single=" << to_string_u128(single)
                          << ", multi=" << to_string_u128(multi) << ")\n";
                return false;
            }
        }
    }

    if (kTargetDiv6 != 75ULL) {
        std::cerr << "Validation failed: target setup mismatch.\n";
        return false;
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) return 1;

    if (options.run_checks) {
        if (!run_validations(options)) return 1;
    }

    const u128 n_limit = n_limit_from_l_limit(options.limit_l);
    const u128 answer = solve_for_n_limit(n_limit,
                                          options.allow_multithreading,
                                          options.requested_threads);
    std::cout << to_string_u128(answer) << '\n';
    return 0;
}
