#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultK = 1'000'000ULL;
constexpr i64 kDefaultX = 1'000'000'000LL;
constexpr u64 kExpectedAnswerForDefaultInput = 141'630'459'461'893'728ULL;

struct Checkpoint {
    u64 k = 0ULL;
    i64 x = 0LL;
    u64 expected = 0ULL;
};

constexpr Checkpoint kCheckpoints[] = {
    {1ULL, 10LL, 41ULL},
    {10ULL, 100LL, 12'492ULL},
    {1'000ULL, 10'000LL, 329'864'366ULL},
};

struct Options {
    u64 k = kDefaultK;
    i64 x = kDefaultX;
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

bool parse_i64_after_prefix(const std::string& arg, const char* prefix, i64& value) {
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

    if (parsed > static_cast<u64>(std::numeric_limits<i64>::max())) {
        return false;
    }

    value = static_cast<i64>(parsed);
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
        if (parse_u64_after_prefix(arg, "--k=", parsed_u64)) {
            options.k = parsed_u64;
            continue;
        }

        i64 parsed_i64 = 0LL;
        if (parse_i64_after_prefix(arg, "--x=", parsed_i64)) {
            options.x = parsed_i64;
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

    if (options.k == 0ULL) {
        std::cerr << "--k must be positive.\n";
        return false;
    }
    if (options.x <= 0LL) {
        std::cerr << "--x must be positive.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const std::size_t workload_units) {
    constexpr std::size_t kMinUnitsForParallel = 1'000ULL;

    if (!allow_multithreading || workload_units < 2ULL || workload_units < kMinUnitsForParallel) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload_units)));
}

i64 floor_div(const i64 a, const i64 b) {
    if (a >= 0LL) {
        return a / b;
    }
    return -(((-a) + b - 1LL) / b);
}

u64 count_interval(const i64 lo, const i64 hi) {
    if (lo > hi) {
        return 0ULL;
    }
    return static_cast<u64>(hi - lo + 1LL);
}

u64 count_angle_a_for_pair(const i64 u, const i64 v, const i64 x_bound) {
    const i64 lo = std::max(-x_bound, v - x_bound);
    const i64 hi = std::min({x_bound, u + x_bound, floor_div(u - 1LL, 2LL)});
    return count_interval(lo, hi);
}

u64 count_angle_b_for_pair(const i64 u, const i64 v, const i64 x_bound) {
    const i64 lo = std::max({-x_bound, v - x_bound, floor_div(u, 2LL) + 1LL});
    const i64 hi = std::min({x_bound, u + x_bound, floor_div(v - 1LL, 2LL)});
    return count_interval(lo, hi);
}

u64 count_angle_c_for_pair(const i64 u, const i64 v, const i64 x_bound) {
    const i64 lo = std::max({-x_bound, v - x_bound, floor_div(v, 2LL) + 1LL});
    const i64 hi = std::min(x_bound, u + x_bound);
    return count_interval(lo, hi);
}

bool valid_triangle_from_sums(const i64 s1, const i64 s2, const i64 s3, const i64 x_bound) {
    if (!(s1 < s2 && s2 < s3)) {
        return false;
    }

    const i64 parity_sum = s1 + s2 + s3;
    if ((parity_sum & 1LL) != 0LL) {
        return false;
    }

    const i64 a = (s1 + s2 - s3) / 2LL;
    const i64 b = (s1 + s3 - s2) / 2LL;
    const i64 c = (s2 + s3 - s1) / 2LL;

    if (a < -x_bound || a > x_bound) {
        return false;
    }
    if (b < -x_bound || b > x_bound) {
        return false;
    }
    if (c < -x_bound || c > x_bound) {
        return false;
    }

    return a < b && b < c;
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

std::vector<int> build_spf(const u64 limit) {
    std::vector<int> spf(static_cast<std::size_t>(limit) + 1ULL, 0);
    if (limit >= 1ULL) {
        spf[1] = 1;
    }

    for (u64 i = 2ULL; i <= limit; ++i) {
        if (spf[static_cast<std::size_t>(i)] == 0) {
            spf[static_cast<std::size_t>(i)] = static_cast<int>(i);
            if (i <= limit / i) {
                for (u64 j = i * i; j <= limit; j += i) {
                    if (spf[static_cast<std::size_t>(j)] == 0) {
                        spf[static_cast<std::size_t>(j)] = static_cast<int>(i);
                    }
                }
            }
        }
    }

    return spf;
}

void collect_divisors_recursive(const std::vector<std::pair<u64, int>>& factors,
                                const std::size_t index,
                                const u64 current,
                                std::vector<u64>& divisors) {
    if (index == factors.size()) {
        divisors.push_back(current);
        return;
    }

    const u64 prime = factors[index].first;
    const int exponent = factors[index].second;

    u64 value = current;
    for (int e = 0; e <= exponent; ++e) {
        collect_divisors_recursive(factors, index + 1ULL, value, divisors);
        if (e != exponent) {
            value *= prime;
        }
    }
}

void factor_2k_squared(const u64 k,
                       const std::vector<int>& spf,
                       std::vector<std::pair<u64, int>>& factors) {
    factors.clear();

    u64 n = k;
    int exponent_two_in_k = 0;
    while ((n & 1ULL) == 0ULL) {
        ++exponent_two_in_k;
        n >>= 1U;
    }

    factors.emplace_back(2ULL, 2 * exponent_two_in_k + 1);

    while (n > 1ULL) {
        const int p = spf[static_cast<std::size_t>(n)];
        int exp_in_k = 0;
        while (n % static_cast<u64>(p) == 0ULL) {
            n /= static_cast<u64>(p);
            ++exp_in_k;
        }
        factors.emplace_back(static_cast<u64>(p), 2 * exp_in_k);
    }
}

void build_pairs_for_k(const u64 k,
                       const i64 x_bound,
                       const std::vector<u64>& divisors,
                       std::vector<std::pair<i64, i64>>& plus_pairs,
                       std::vector<std::pair<i64, i64>>& minus_pairs,
                       std::vector<std::pair<i64, i64>>& plus_out,
                       std::vector<std::pair<i64, i64>>& plus_in,
                       std::vector<std::pair<i64, i64>>& minus_out,
                       std::vector<std::pair<i64, i64>>& minus_in) {
    plus_pairs.clear();
    minus_pairs.clear();
    plus_out.clear();
    plus_in.clear();
    minus_out.clear();
    minus_in.clear();

    const i64 kk = static_cast<i64>(k);
    const i64 n = static_cast<i64>(2ULL * k * k);
    const i64 sum_bound = 2LL * x_bound;

    plus_pairs.reserve(divisors.size() * 2ULL);
    minus_pairs.reserve(divisors.size() * 2ULL);

    for (const u64 d_u64 : divisors) {
        const i64 d = static_cast<i64>(d_u64);
        for (const int sign : {-1, 1}) {
            const i64 p = static_cast<i64>(sign) * d;
            const i64 q = -n / p;

            const i64 u_plus = p + kk;
            const i64 v_plus = q - kk;
            if (u_plus < v_plus &&
                -sum_bound <= u_plus && u_plus <= sum_bound &&
                -sum_bound <= v_plus && v_plus <= sum_bound) {
                plus_pairs.emplace_back(u_plus, v_plus);
            }

            const i64 u_minus = p - kk;
            const i64 v_minus = q + kk;
            if (u_minus < v_minus &&
                -sum_bound <= u_minus && u_minus <= sum_bound &&
                -sum_bound <= v_minus && v_minus <= sum_bound) {
                minus_pairs.emplace_back(u_minus, v_minus);
            }
        }
    }

    std::sort(plus_pairs.begin(), plus_pairs.end());
    plus_pairs.erase(std::unique(plus_pairs.begin(), plus_pairs.end()), plus_pairs.end());

    std::sort(minus_pairs.begin(), minus_pairs.end());
    minus_pairs.erase(std::unique(minus_pairs.begin(), minus_pairs.end()), minus_pairs.end());

    plus_out = plus_pairs;
    minus_out = minus_pairs;

    plus_in.reserve(plus_pairs.size());
    for (const auto& [u, v] : plus_pairs) {
        plus_in.emplace_back(v, u);
    }
    std::sort(plus_in.begin(), plus_in.end());

    minus_in.reserve(minus_pairs.size());
    for (const auto& [u, v] : minus_pairs) {
        minus_in.emplace_back(v, u);
    }
    std::sort(minus_in.begin(), minus_in.end());
}

u128 solve_range(const u64 thread_index,
                 const u64 thread_count,
                 const u64 k_limit,
                 const i64 x_bound,
                 const std::vector<int>& spf) {
    u128 subtotal = 0U;

    std::vector<std::pair<u64, int>> factors;
    std::vector<u64> divisors;

    std::vector<std::pair<i64, i64>> plus_pairs;
    std::vector<std::pair<i64, i64>> minus_pairs;
    std::vector<std::pair<i64, i64>> plus_out;
    std::vector<std::pair<i64, i64>> plus_in;
    std::vector<std::pair<i64, i64>> minus_out;
    std::vector<std::pair<i64, i64>> minus_in;

    factors.reserve(16ULL);
    divisors.reserve(2048ULL);

    for (u64 k = 1ULL + thread_index; k <= k_limit; k += thread_count) {
        factor_2k_squared(k, spf, factors);

        divisors.clear();
        collect_divisors_recursive(factors, 0ULL, 1ULL, divisors);

        build_pairs_for_k(k,
                          x_bound,
                          divisors,
                          plus_pairs,
                          minus_pairs,
                          plus_out,
                          plus_in,
                          minus_out,
                          minus_in);

        u128 count_a = 0U;
        u128 count_b = 0U;
        u128 count_c = 0U;

        for (const auto& [u, v] : plus_out) {
            count_a += static_cast<u128>(count_angle_a_for_pair(u, v, x_bound));
            count_c += static_cast<u128>(count_angle_c_for_pair(u, v, x_bound));
        }

        for (const auto& [u, v] : minus_out) {
            count_b += static_cast<u128>(count_angle_b_for_pair(u, v, x_bound));
        }

        u128 overlap_ab = 0U;
        {
            std::size_t i = 0ULL;
            std::size_t j = 0ULL;
            while (i < plus_out.size() && j < minus_out.size()) {
                const i64 key_plus = plus_out[i].first;
                const i64 key_minus = minus_out[j].first;

                if (key_plus < key_minus) {
                    ++i;
                } else if (key_minus < key_plus) {
                    ++j;
                } else {
                    const i64 x = key_plus;
                    const i64 y = plus_out[i].second;
                    const i64 z = minus_out[j].second;
                    if (valid_triangle_from_sums(x, y, z, x_bound)) {
                        ++overlap_ab;
                    }
                    ++i;
                    ++j;
                }
            }
        }

        u128 overlap_ac = 0U;
        {
            std::size_t i = 0ULL;
            std::size_t j = 0ULL;
            while (i < plus_in.size() && j < plus_out.size()) {
                const i64 key_in = plus_in[i].first;
                const i64 key_out = plus_out[j].first;

                if (key_in < key_out) {
                    ++i;
                } else if (key_out < key_in) {
                    ++j;
                } else {
                    const i64 y = key_in;
                    const i64 x = plus_in[i].second;
                    const i64 z = plus_out[j].second;
                    if (valid_triangle_from_sums(x, y, z, x_bound)) {
                        ++overlap_ac;
                    }
                    ++i;
                    ++j;
                }
            }
        }

        u128 overlap_bc = 0U;
        {
            std::size_t i = 0ULL;
            std::size_t j = 0ULL;
            while (i < plus_in.size() && j < minus_in.size()) {
                const i64 key_plus = plus_in[i].first;
                const i64 key_minus = minus_in[j].first;

                if (key_plus < key_minus) {
                    ++i;
                } else if (key_minus < key_plus) {
                    ++j;
                } else {
                    const i64 z = key_plus;
                    const i64 y = plus_in[i].second;
                    const i64 x = minus_in[j].second;
                    if (valid_triangle_from_sums(x, y, z, x_bound)) {
                        ++overlap_bc;
                    }
                    ++i;
                    ++j;
                }
            }
        }

        subtotal += count_a + count_b + count_c - overlap_ab - overlap_ac - overlap_bc;
    }

    return subtotal;
}

u128 solve_f(const u64 k_limit,
             const i64 x_bound,
             const bool allow_multithreading,
             const unsigned requested_threads,
             const std::vector<int>& spf) {
    const unsigned thread_count =
        choose_thread_count(allow_multithreading, requested_threads, static_cast<std::size_t>(k_limit));

    if (thread_count == 1U) {
        return solve_range(0ULL, 1ULL, k_limit, x_bound, spf);
    }

    std::vector<u128> partial(thread_count, 0U);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    for (unsigned t = 0U; t < thread_count; ++t) {
        workers.emplace_back([&, t]() {
            partial[t] = solve_range(static_cast<u64>(t),
                                     static_cast<u64>(thread_count),
                                     k_limit,
                                     x_bound,
                                     spf);
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    u128 total = 0U;
    for (const u128 value : partial) {
        total += value;
    }
    return total;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        std::cerr << "Usage: " << argv[0]
                  << " [--k=<value>] [--x=<value>] [--threads=<value>]"
                  << " [--single-thread] [--skip-checkpoints]\n";
        return 1;
    }

    u64 max_k_needed = options.k;
    if (options.run_checkpoints) {
        for (const Checkpoint& cp : kCheckpoints) {
            max_k_needed = std::max(max_k_needed, cp.k);
        }
    }

    const std::vector<int> spf = build_spf(max_k_needed);

    if (options.run_checkpoints) {
        for (const Checkpoint& cp : kCheckpoints) {
            const u128 actual =
                solve_f(cp.k, cp.x, options.allow_multithreading, options.requested_threads, spf);
            const u128 expected = static_cast<u128>(cp.expected);
            std::cout << "Checkpoint F(" << cp.k << ", " << cp.x << ") = "
                      << to_string_u128(actual) << " (expected " << cp.expected << ")\n";
            if (actual != expected) {
                std::cerr << "Checkpoint mismatch.\n";
                return 1;
            }
        }
    }

    const u128 answer =
        solve_f(options.k, options.x, options.allow_multithreading, options.requested_threads, spf);

    std::cout << "F(" << options.k << ", " << options.x << ") = " << to_string_u128(answer)
              << '\n';
    if (options.k == kDefaultK && options.x == kDefaultX) {
        const u128 expected = static_cast<u128>(kExpectedAnswerForDefaultInput);
        if (answer != expected) {
            std::cerr << "Default-case mismatch: expected "
                      << kExpectedAnswerForDefaultInput << ".\n";
            return 1;
        }
    }

    return 0;
}
