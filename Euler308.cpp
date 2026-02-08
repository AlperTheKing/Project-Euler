#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = __uint128_t;

struct Options {
    u64 target_index = 10001ULL;
    bool run_checkpoints = true;
};

bool parse_u64_after_prefix(const std::string& arg,
                            const std::string& prefix,
                            u64& value) {
    if (arg.rfind(prefix, 0U) != 0U) {
        return false;
    }

    const std::string tail = arg.substr(prefix.size());
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
            throw std::overflow_error("Argument overflow");
        }
        parsed = parsed * 10ULL + digit;
    }

    value = parsed;
    return true;
}

bool parse_arguments(const int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--target=", options.target_index)) {
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

std::string to_string_u128(u128 value) {
    if (value == 0U) {
        return "0";
    }

    std::string digits;
    while (value > 0U) {
        const int digit = static_cast<int>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }

    std::reverse(digits.begin(), digits.end());
    return digits;
}

u64 estimate_upper_bound_for_nth_prime(const u64 n) {
    if (n < 6ULL) {
        return 15ULL;
    }

    const long double nd = static_cast<long double>(n);
    const long double estimate =
        nd * (std::log(nd) + std::log(std::log(nd))) + 32.0L;

    if (!std::isfinite(static_cast<double>(estimate)) || estimate <= 0.0L) {
        throw std::runtime_error("Prime bound estimate failed");
    }

    const u64 bound = static_cast<u64>(std::ceil(estimate));
    return std::max<u64>(bound, 15ULL);
}

struct SieveData {
    std::vector<u32> primes;
    std::vector<u32> spf;
};

SieveData build_linear_sieve(const u64 limit) {
    if (limit > static_cast<u64>(std::numeric_limits<u32>::max())) {
        throw std::runtime_error("Limit too large for 32-bit SPF table");
    }

    SieveData out;
    out.spf.assign(static_cast<std::size_t>(limit) + 1ULL, 0U);
    out.primes.reserve(static_cast<std::size_t>(limit / 10ULL + 16ULL));

    for (u32 i = 2U; i <= limit; ++i) {
        if (out.spf[static_cast<std::size_t>(i)] == 0U) {
            out.spf[static_cast<std::size_t>(i)] = i;
            out.primes.push_back(i);
        }

        for (const u32 p : out.primes) {
            const u64 composite = static_cast<u64>(p) * static_cast<u64>(i);
            if (composite > limit) {
                break;
            }

            out.spf[static_cast<std::size_t>(composite)] = p;
            if (p == out.spf[static_cast<std::size_t>(i)]) {
                break;
            }
        }
    }

    return out;
}

SieveData build_sieve_with_prime_count(const u64 target_index) {
    u64 limit = estimate_upper_bound_for_nth_prime(target_index);

    while (true) {
        SieveData sieve = build_linear_sieve(limit);
        if (sieve.primes.size() >= target_index) {
            return sieve;
        }
        if (limit > std::numeric_limits<u64>::max() / 2ULL) {
            throw std::runtime_error("Prime sieve limit overflow");
        }
        limit *= 2ULL;
    }
}

u128 floor_sum_upto(const u64 n, u64 r) {
    if (r == 0ULL) {
        return 0U;
    }
    if (r > n) {
        r = n;
    }

    u128 sum = 0U;
    u64 left = 1ULL;
    while (left <= r) {
        const u64 q = n / left;
        u64 right = n / q;
        if (right > r) {
            right = r;
        }

        sum += static_cast<u128>(q) * static_cast<u128>(right - left + 1ULL);
        left = right + 1ULL;
    }

    return sum;
}

u128 floor_sum_range(const u64 n, const u64 l, const u64 r) {
    if (l > r || r == 0ULL) {
        return 0U;
    }
    return floor_sum_upto(n, r) - floor_sum_upto(n, l - 1ULL);
}

u128 candidate_delta_to_next_start(const u64 n, const u64 spf_n) {
    const u64 max_proper_divisor = n / spf_n;

    const u64 non_divisor_count =
        (n > max_proper_divisor + 1ULL) ? (n - max_proper_divisor - 1ULL) : 0ULL;
    const u128 non_divisor_floor_sum =
        floor_sum_range(n, max_proper_divisor + 1ULL, n - 1ULL);

    // Candidate processing decomposition:
    // setup in S11, all non-divisor tests, and terminal divisor branch.
    const u128 setup = static_cast<u128>(2ULL * n - 1ULL);
    const u128 non_divisor_total =
        static_cast<u128>(6ULL * n + 2ULL) * static_cast<u128>(non_divisor_count) +
        2U * non_divisor_floor_sum;
    const u128 terminal =
        static_cast<u128>(5ULL * n + max_proper_divisor + 2ULL * (n / max_proper_divisor) +
                          2ULL);

    return setup + non_divisor_total + terminal;
}

u128 prime_hit_offset_from_start(const u64 prime_candidate) {
    const u128 non_divisor_floor_sum = floor_sum_range(prime_candidate, 2ULL, prime_candidate - 1ULL);

    const u128 setup = static_cast<u128>(2ULL * prime_candidate - 1ULL);
    const u128 non_divisor_total =
        static_cast<u128>(6ULL * prime_candidate + 2ULL) *
            static_cast<u128>(prime_candidate - 2ULL) +
        2U * non_divisor_floor_sum;
    const u128 prime_branch_to_hit = static_cast<u128>(6ULL * prime_candidate + 2ULL);

    return setup + non_divisor_total + prime_branch_to_hit;
}

u128 solve_iterations_to_prime_index(const u64 target_index) {
    if (target_index == 0ULL) {
        throw std::invalid_argument("--target must be positive");
    }

    const SieveData sieve = build_sieve_with_prime_count(target_index);
    const u64 prime_target =
        static_cast<u64>(sieve.primes[static_cast<std::size_t>(target_index - 1ULL)]);

    // T_n is the first step index where the machine is in canonical start-state
    // for candidate n. T_2 = 2 from seed 2 after applying rules 12 and 14 once each.
    u128 start_step = 2U;
    u64 found_primes = 0ULL;

    for (u64 n = 2ULL; n <= prime_target; ++n) {
        const u64 spf_n = static_cast<u64>(sieve.spf[static_cast<std::size_t>(n)]);
        const bool is_prime = (spf_n == n);

        if (is_prime) {
            const u128 hit_step = start_step + prime_hit_offset_from_start(n);
            ++found_primes;
            if (found_primes == target_index) {
                return hit_step;
            }
        }

        start_step += candidate_delta_to_next_start(n, spf_n);
    }

    throw std::runtime_error("Prime target was not reached");
}

bool is_prime_small(const u64 x) {
    if (x < 2ULL) {
        return false;
    }
    if ((x % 2ULL) == 0ULL) {
        return x == 2ULL;
    }
    for (u64 d = 3ULL; d * d <= x; d += 2ULL) {
        if ((x % d) == 0ULL) {
            return false;
        }
    }
    return true;
}

struct FractranState {
    u64 a = 1ULL;  // exponent of 2
    u64 b = 0ULL;  // exponent of 3
    u64 c = 0ULL;  // exponent of 5
    u64 d = 0ULL;  // exponent of 7
    u64 e = 0ULL;  // exponent of 11
    u64 f = 0ULL;  // exponent of 13
    u64 g = 0ULL;  // exponent of 17
    u64 h = 0ULL;  // exponent of 19
    u64 i = 0ULL;  // exponent of 23
    u64 j = 0ULL;  // exponent of 29
};

bool is_pure_power_of_two(const FractranState& s) {
    return s.b == 0ULL && s.c == 0ULL && s.d == 0ULL && s.e == 0ULL && s.f == 0ULL &&
           s.g == 0ULL && s.h == 0ULL && s.i == 0ULL && s.j == 0ULL;
}

void apply_one_fractran_step(FractranState& s) {
    if (s.d > 0ULL && s.f > 0ULL) {
        --s.d;
        --s.f;
        ++s.g;
        return;
    }
    if (s.c > 0ULL && s.g > 0ULL) {
        ++s.a;
        ++s.b;
        --s.c;
        ++s.f;
        --s.g;
        return;
    }
    if (s.b > 0ULL && s.g > 0ULL) {
        --s.b;
        --s.g;
        ++s.h;
        return;
    }
    if (s.a > 0ULL && s.h > 0ULL) {
        --s.a;
        --s.h;
        ++s.i;
        return;
    }
    if (s.b > 0ULL && s.e > 0ULL) {
        --s.b;
        --s.e;
        ++s.j;
        return;
    }
    if (s.j > 0ULL) {
        ++s.d;
        ++s.e;
        --s.j;
        return;
    }
    if (s.i > 0ULL) {
        ++s.c;
        ++s.h;
        --s.i;
        return;
    }
    if (s.h > 0ULL) {
        ++s.d;
        ++s.e;
        --s.h;
        return;
    }
    if (s.g > 0ULL) {
        --s.g;
        return;
    }
    if (s.f > 0ULL) {
        ++s.e;
        --s.f;
        return;
    }
    if (s.e > 0ULL) {
        --s.e;
        ++s.f;
        return;
    }
    if (s.a > 0ULL) {
        --s.a;
        ++s.b;
        ++s.c;
        return;
    }
    if (s.d > 0ULL) {
        --s.d;
        return;
    }

    ++s.c;
    ++s.e;
}

std::vector<std::pair<u64, u64>> brute_prime_hits(const std::size_t count) {
    std::vector<std::pair<u64, u64>> hits;
    hits.reserve(count);

    FractranState state;
    u64 steps = 0ULL;

    while (hits.size() < count) {
        if (is_pure_power_of_two(state) && is_prime_small(state.a)) {
            hits.push_back({state.a, steps});
        }

        apply_one_fractran_step(state);
        ++steps;
        if (steps > 2000000000ULL) {
            throw std::runtime_error("Brute-force checkpoint exceeded step cap");
        }
    }

    return hits;
}

bool check_equal_u64(const u64 actual, const u64 expected, const std::string& label) {
    if (actual == expected) {
        return true;
    }
    std::cerr << "Checkpoint failed for " << label << ": expected " << expected
              << ", got " << actual << '\n';
    return false;
}

bool check_equal_u128(const u128 actual,
                      const u128 expected,
                      const std::string& label) {
    if (actual == expected) {
        return true;
    }
    std::cerr << "Checkpoint failed for " << label << ": expected "
              << to_string_u128(expected) << ", got " << to_string_u128(actual)
              << '\n';
    return false;
}

bool run_checkpoints() {
    const std::vector<std::pair<u64, u64>> brute10 = brute_prime_hits(10U);
    const std::vector<std::pair<u64, u64>> expected10 = {
        {2ULL, 19ULL},   {3ULL, 69ULL},    {5ULL, 281ULL},   {7ULL, 710ULL},
        {11ULL, 2375ULL}, {13ULL, 3893ULL}, {17ULL, 8102ULL}, {19ULL, 11361ULL},
        {23ULL, 19268ULL}, {29ULL, 36981ULL}};

    if (brute10.size() != expected10.size()) {
        std::cerr << "Checkpoint failed: brute sequence size mismatch\n";
        return false;
    }
    for (std::size_t i = 0; i < expected10.size(); ++i) {
        if (!check_equal_u64(brute10[i].first, expected10[i].first,
                             "prime exponent #" + std::to_string(i + 1U))) {
            return false;
        }
        if (!check_equal_u64(brute10[i].second, expected10[i].second,
                             "iteration count #" + std::to_string(i + 1U))) {
            return false;
        }
    }

    const std::vector<std::pair<u64, u64>> brute30 = brute_prime_hits(30U);
    const std::vector<u64> targets = {1ULL, 2ULL, 3ULL, 10ULL, 20ULL, 30ULL};
    for (const u64 target : targets) {
        const u128 formula = solve_iterations_to_prime_index(target);
        const u128 expected =
            static_cast<u128>(brute30[static_cast<std::size_t>(target - 1ULL)].second);
        if (!check_equal_u128(formula, expected,
                              "formula-vs-brute for target=" + std::to_string(target))) {
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        Options options;
        if (!parse_arguments(argc, argv, options)) {
            return 1;
        }
        if (options.target_index == 0ULL) {
            std::cerr << "--target must be positive.\n";
            return 1;
        }

        if (options.run_checkpoints && !run_checkpoints()) {
            return 1;
        }

        const u128 answer = solve_iterations_to_prime_index(options.target_index);
        std::cout << to_string_u128(answer) << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
