#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace {

using boost::multiprecision::cpp_int;
using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr i64 kMod = 1'000'000'007LL;
constexpr u64 kDefaultK = 10'000ULL;

struct Options {
    u64 k = kDefaultK;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

u64 cube_u64(const u64 n) {
    const u128 x = static_cast<u128>(n);
    return static_cast<u64>(x * x * x);
}

u64 fourth_u64(const u64 n) {
    const u128 x = static_cast<u128>(n);
    return static_cast<u64>(x * x * x * x);
}

cpp_int binom_small_k(u64 n, u64 k) {
    if (k > n) {
        return 0;
    }
    if (k > n - k) {
        k = n - k;
    }

    cpp_int out = 1;
    for (u64 i = 1; i <= k; ++i) {
        out *= (n - k + i);
        out /= i;
    }
    return out;
}

cpp_int count_0_to_9(const u64 digits, const u64 sum) {
    if (sum > 9ULL * digits) {
        return 0;
    }
    if (digits == 0ULL) {
        return (sum == 0ULL) ? cpp_int(1) : cpp_int(0);
    }

    cpp_int total = 0;
    const u64 max_j = sum / 10ULL;
    for (u64 j = 0; j <= max_j; ++j) {
        const u64 remain = sum - 10ULL * j;

        cpp_int term = binom_small_k(digits, j);
        term *= binom_small_k(digits + remain - 1ULL, remain);

        if ((j & 1ULL) == 0ULL) {
            total += term;
        } else {
            total -= term;
        }
    }

    return total;
}

cpp_int count_len_with_sum(const u64 len, const u64 sum) {
    if (len == 0ULL || sum == 0ULL) {
        return 0;
    }
    if (sum > 9ULL * len) {
        return 0;
    }

    const u64 deficit = 9ULL * len - sum;
    const u64 max_first_deficit = std::min<u64>(8ULL, deficit);

    cpp_int out = 0;
    for (u64 b0 = 0ULL; b0 <= max_first_deficit; ++b0) {
        out += count_0_to_9(len - 1ULL, deficit - b0);
    }
    return out;
}

i64 mod_pow10(u64 exp) {
    i64 base = 10LL;
    i64 out = 1LL;

    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) {
            out = static_cast<i64>((static_cast<__int128>(out) * base) % kMod);
        }
        base = static_cast<i64>((static_cast<__int128>(base) * base) % kMod);
        exp >>= 1ULL;
    }

    return out;
}

void append_digit_mod(i64& value_mod, const int digit) {
    value_mod = static_cast<i64>((static_cast<__int128>(value_mod) * 10LL + digit) % kMod);
}

void append_nines_mod(i64& value_mod, const u64 count) {
    if (count == 0ULL) {
        return;
    }

    const i64 p10 = mod_pow10(count);
    // 9 repeated count times = 10^count - 1.
    value_mod = static_cast<i64>(
        (static_cast<__int128>(value_mod) * p10 + p10 - 1LL + kMod) % kMod);
}

struct RankedNumber {
    u64 len = 0ULL;
    cpp_int rank_inside_len = 0;
};

RankedNumber locate_length_and_rank(const u64 digit_sum, const u64 occurrence_rank) {
    const cpp_int target = occurrence_rank;

    u64 len = (digit_sum + 8ULL) / 9ULL;
    cpp_int prefix_count = 0;

    while (true) {
        const cpp_int cnt_here = count_len_with_sum(len, digit_sum);
        if (prefix_count + cnt_here >= target) {
            return RankedNumber{len, target - prefix_count};
        }
        prefix_count += cnt_here;
        ++len;
    }
}

i64 unrank_mod(const u64 len, const u64 digit_sum, cpp_int rank_inside_len) {
    const u64 initial_deficit = 9ULL * len - digit_sum;

    i64 value_mod = 0LL;
    u64 remaining_digits = len;
    u64 remaining_deficit = initial_deficit;

    {
        const u64 max_deficit_here = std::min<u64>(8ULL, remaining_deficit);
        bool chosen = false;

        for (int deficit_digit = static_cast<int>(max_deficit_here); deficit_digit >= 0;
             --deficit_digit) {
            const cpp_int cnt =
                count_0_to_9(remaining_digits - 1ULL, remaining_deficit - deficit_digit);

            if (rank_inside_len > cnt) {
                rank_inside_len -= cnt;
                continue;
            }

            append_digit_mod(value_mod, 9 - deficit_digit);
            remaining_deficit -= static_cast<u64>(deficit_digit);
            --remaining_digits;
            chosen = true;
            break;
        }

        if (!chosen) {
            throw std::runtime_error("Failed to choose the first digit during unranking.");
        }
    }

    while (remaining_digits > 0ULL) {
        if (remaining_deficit == 0ULL) {
            append_nines_mod(value_mod, remaining_digits);
            break;
        }

        const cpp_int total_here = count_0_to_9(remaining_digits, remaining_deficit);
        const cpp_int zero_branch = count_0_to_9(remaining_digits - 1ULL, remaining_deficit);
        const cpp_int nonzero_prefix = total_here - zero_branch;

        if (rank_inside_len > nonzero_prefix) {
            u64 low = 1ULL;
            u64 high = remaining_digits;
            u64 best = 1ULL;

            while (low <= high) {
                const u64 mid = low + ((high - low) >> 1ULL);
                const cpp_int suffix = count_0_to_9(remaining_digits - mid, remaining_deficit);
                const cpp_int removed = total_here - suffix;

                if (rank_inside_len > removed) {
                    best = mid;
                    low = mid + 1ULL;
                } else {
                    if (mid == 0ULL) {
                        break;
                    }
                    high = mid - 1ULL;
                }
            }

            const cpp_int suffix = count_0_to_9(remaining_digits - best, remaining_deficit);
            const cpp_int removed = total_here - suffix;
            rank_inside_len -= removed;

            append_nines_mod(value_mod, best);
            remaining_digits -= best;
            continue;
        }

        bool chosen = false;
        const u64 max_deficit_here = std::min<u64>(9ULL, remaining_deficit);
        for (int deficit_digit = static_cast<int>(max_deficit_here); deficit_digit >= 0;
             --deficit_digit) {
            const cpp_int cnt =
                count_0_to_9(remaining_digits - 1ULL, remaining_deficit - deficit_digit);

            if (rank_inside_len > cnt) {
                rank_inside_len -= cnt;
                continue;
            }

            append_digit_mod(value_mod, 9 - deficit_digit);
            remaining_deficit -= static_cast<u64>(deficit_digit);
            --remaining_digits;
            chosen = true;
            break;
        }

        if (!chosen) {
            throw std::runtime_error("Failed to choose a digit during unranking.");
        }
    }

    return value_mod;
}

i64 f_mod(const u64 digit_sum, const u64 occurrence_rank) {
    const RankedNumber located = locate_length_and_rank(digit_sum, occurrence_rank);
    return unrank_mod(located.len, digit_sum, located.rank_inside_len);
}

int digit_sum_small(u64 x) {
    int out = 0;
    while (x > 0ULL) {
        out += static_cast<int>(x % 10ULL);
        x /= 10ULL;
    }
    return out;
}

u64 brute_f_small(const int target_sum, const int target_occurrence) {
    int seen = 0;
    for (u64 x = 1ULL;; ++x) {
        if (digit_sum_small(x) != target_sum) {
            continue;
        }
        ++seen;
        if (seen == target_occurrence) {
            return x;
        }
    }
}

unsigned choose_thread_count(const bool allow_multithreading,
                             const unsigned requested_threads,
                             const u64 k) {
    if (!allow_multithreading || k <= 1ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(k)));
}

i64 solve_S_mod(const u64 k, const unsigned threads) {
    if (k == 0ULL) {
        return 0LL;
    }

    if (threads <= 1U) {
        i64 total = 0LL;
        for (u64 n = 1ULL; n <= k; ++n) {
            const u64 s = cube_u64(n);
            const u64 m = fourth_u64(n);
            total += f_mod(s, m);
            if (total >= kMod) {
                total -= kMod;
            }
        }
        return total;
    }

    std::vector<i64> partial(threads, 0LL);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    for (unsigned tid = 0U; tid < threads; ++tid) {
        const u64 begin = 1ULL + (k * tid) / threads;
        const u64 end = (k * (tid + 1U)) / threads;

        pool.emplace_back([&, tid, begin, end]() {
            i64 local = 0LL;
            for (u64 n = begin; n <= end; ++n) {
                const u64 s = cube_u64(n);
                const u64 m = fourth_u64(n);
                local += f_mod(s, m);
                if (local >= kMod) {
                    local -= kMod;
                }
            }
            partial[tid] = local;
        });
    }

    for (std::thread& th : pool) {
        th.join();
    }

    i64 total = 0LL;
    for (const i64 part : partial) {
        total += part;
        total %= kMod;
    }

    return total;
}

bool parse_u64(const std::string& s, u64& out) {
    if (s.empty()) {
        return false;
    }

    u64 value = 0ULL;
    for (const char c : s) {
        if (c < '0' || c > '9') {
            return false;
        }
        const u64 digit = static_cast<u64>(c - '0');
        if (value > (std::numeric_limits<u64>::max() - digit) / 10ULL) {
            return false;
        }
        value = value * 10ULL + digit;
    }

    out = value;
    return true;
}

bool parse_unsigned(const std::string& s, unsigned& out) {
    u64 parsed = 0ULL;
    if (!parse_u64(s, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }
    out = static_cast<unsigned>(parsed);
    return true;
}

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& out) {
    const std::string p(prefix);
    if (arg.rfind(p, 0U) != 0U) {
        return false;
    }
    return parse_u64(arg.substr(p.size()), out);
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& out) {
    const std::string p(prefix);
    if (arg.rfind(p, 0U) != 0U) {
        return false;
    }
    return parse_unsigned(arg.substr(p.size()), out);
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

        u64 parsed_k = 0ULL;
        if (parse_u64_after_prefix(arg, "--k=", parsed_k)) {
            options.k = parsed_k;
            continue;
        }

        unsigned parsed_threads = 0U;
        if (parse_unsigned_after_prefix(arg, "--threads=", parsed_threads)) {
            options.requested_threads = parsed_threads;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    return true;
}

bool run_checkpoints(const unsigned threads) {
    struct FCheckpoint {
        u64 sum;
        u64 occ;
        u64 expected;
    };

    const std::vector<FCheckpoint> f_checks = {
        {10ULL, 1ULL, 19ULL},
        {10ULL, 10ULL, 109ULL},
        {10ULL, 100ULL, 1423ULL},
    };

    for (const FCheckpoint cp : f_checks) {
        const u64 brute = brute_f_small(static_cast<int>(cp.sum), static_cast<int>(cp.occ));
        if (brute != cp.expected) {
            std::cerr << "Internal brute checkpoint mismatch for f(" << cp.sum << ',' << cp.occ
                      << "): expected " << cp.expected << ", brute got " << brute << "\n";
            return false;
        }

        const i64 got_mod = f_mod(cp.sum, cp.occ);
        if (got_mod != static_cast<i64>(cp.expected % kMod)) {
            std::cerr << "f checkpoint failed for f(" << cp.sum << ',' << cp.occ << "): expected "
                      << cp.expected << ", got mod value " << got_mod << "\n";
            return false;
        }
    }

    const i64 s3 = solve_S_mod(3ULL, 1U);
    if (s3 != 7128LL) {
        std::cerr << "S(3) checkpoint failed: expected 7128, got " << s3 << "\n";
        return false;
    }

    const i64 s10 = solve_S_mod(10ULL, 1U);
    if (s10 != 32'287'064LL) {
        std::cerr << "S(10) checkpoint failed: expected 32287064, got " << s10 << "\n";
        return false;
    }

    if (threads > 1U) {
        const i64 single = solve_S_mod(120ULL, 1U);
        const i64 multi = solve_S_mod(120ULL, threads);
        if (single != multi) {
            std::cerr << "Thread consistency checkpoint failed: single=" << single
                      << ", multi=" << multi << "\n";
            return false;
        }
    }

    std::cout << "All checkpoints passed.\n";
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const unsigned threads =
        choose_thread_count(options.allow_multithreading, options.requested_threads, options.k);

    if (options.run_checkpoints && !run_checkpoints(threads)) {
        return 1;
    }

    const i64 answer = solve_S_mod(options.k, threads);
    std::cout << answer << '\n';
    return 0;
}
