#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Options {
    u64 limit = 10'000'000ULL;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0;
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

bool parse_unsigned_after_prefix(const std::string& arg, const std::string& prefix, unsigned& value) {
    u64 parsed = 0ULL;
    if (!parse_u64_after_prefix(arg, prefix, parsed) || parsed == 0ULL ||
        parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
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
        if (arg == "--no-mt") {
            options.allow_multithreading = false;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--limit=", options.limit)) {
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

u64 isqrt_u64(const u64 n) {
    u64 low = 0ULL;
    u64 high = std::min<u64>(n, 1ULL << 32);
    while (low < high) {
        const u64 mid = low + (high - low + 1ULL) / 2ULL;
        if (mid <= n / mid) {
            low = mid;
        } else {
            high = mid - 1ULL;
        }
    }
    return low;
}

u128 arithmetic_series_sum(const u64 count) {
    return static_cast<u128>(count) * static_cast<u128>(count + 1ULL) / 2U;
}

u128 contribution(const u64 limit, const u64 base_sum) {
    const u64 copies = limit / base_sum;
    return static_cast<u128>(base_sum) * arithmetic_series_sum(copies);
}

u128 solve_range(const u64 limit, const u64 n_begin, const u64 n_end) {
    u128 total = 0;

    for (u64 n = n_begin; n < n_end; ++n) {
        const i64 nn = static_cast<i64>(n);

        // Case 1: A = k * 2mn, B = k * (m^2 - n^2), C = k * (m^2 + n^2).
        for (u64 m = n + 1ULL; m <= 5ULL * n; ++m) {
            if (((m - n) & 1ULL) == 0ULL || std::gcd(m, n) != 1ULL) {
                continue;
            }

            const i64 mm = static_cast<i64>(m);
            const i64 b = 5LL * mm * nn - mm * mm - nn * nn;
            const i64 c = mm * mm - 4LL * mm * nn + nn * nn;
            if (b <= 0LL || c <= 0LL) {
                continue;
            }

            const u64 sum_plus = static_cast<u64>(nn * (5LL * mm - nn));
            if (sum_plus <= limit) {
                total += contribution(limit, sum_plus);
            }

            const i64 a_minus = 4LL * mm * nn - mm * mm;
            if (a_minus > 0LL) {
                const u64 sum_minus = static_cast<u64>(mm * (5LL * nn - mm));
                if (sum_minus <= limit) {
                    total += contribution(limit, sum_minus);
                }
            }
        }

        // Case 2: A = 2k * (m^2 - n^2), B = 4k mn, C = 2k * (m^2 + n^2).
        for (u64 m = n + 1ULL; m <= 2ULL * n; ++m) {
            if (((m - n) & 1ULL) == 0ULL || std::gcd(m, n) != 1ULL) {
                continue;
            }

            const i64 mm = static_cast<i64>(m);
            const i64 b = 3LL * mm * mm - 7LL * nn * nn;
            const i64 c = 2LL * (3LL * nn * nn - mm * mm);
            if (b <= 0LL || c <= 0LL) {
                continue;
            }

            const i64 sum_plus = 4LL * mm * mm + 2LL * mm * nn - 6LL * nn * nn;
            if (sum_plus > 0LL && static_cast<u64>(sum_plus) <= limit) {
                total += contribution(limit, static_cast<u64>(sum_plus));
            }

            const i64 a_minus = 3LL * mm * mm - 2LL * mm * nn - 5LL * nn * nn;
            if (a_minus > 0LL) {
                const i64 sum_minus = 4LL * mm * mm - 2LL * mm * nn - 6LL * nn * nn;
                if (sum_minus > 0LL && static_cast<u64>(sum_minus) <= limit) {
                    total += contribution(limit, static_cast<u64>(sum_minus));
                }
            }
        }
    }

    return total;
}

unsigned resolve_thread_count(const Options& options, const u64 max_n) {
    if (!options.allow_multithreading || max_n <= 1ULL) {
        return 1U;
    }
    if (options.requested_threads != 0U) {
        return std::min<unsigned>(options.requested_threads, static_cast<unsigned>(max_n));
    }
    const unsigned detected = std::thread::hardware_concurrency();
    const unsigned fallback = detected == 0U ? 1U : detected;
    return std::min<unsigned>(fallback, static_cast<unsigned>(max_n));
}

u128 solve(const u64 limit, const unsigned thread_count) {
    if (limit == 0ULL) {
        return 0;
    }

    const u64 max_n = isqrt_u64(limit) + 2ULL;
    if (thread_count <= 1U || max_n <= 1ULL) {
        return solve_range(limit, 1ULL, max_n + 1ULL);
    }

    std::vector<u128> partials(thread_count, 0);
    std::vector<std::thread> workers;
    workers.reserve(thread_count);

    const u64 block = (max_n + static_cast<u64>(thread_count) - 1ULL) / static_cast<u64>(thread_count);
    for (unsigned tid = 0; tid < thread_count; ++tid) {
        const u64 begin = 1ULL + static_cast<u64>(tid) * block;
        const u64 end = std::min<u64>(max_n + 1ULL, begin + block);
        workers.emplace_back([&, tid, begin, end]() {
            if (begin < end) {
                partials[tid] = solve_range(limit, begin, end);
            }
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    u128 total = 0;
    for (const u128 part : partials) {
        total += part;
    }
    return total;
}

u128 brute_force(const u64 limit) {
    u128 total = 0;
    for (u64 a = 1ULL; a <= limit; ++a) {
        for (u64 b = 1ULL; a + b <= limit; ++b) {
            for (u64 c = 1ULL; a + b + c <= limit; ++c) {
                const u128 lhs = static_cast<u128>(a) * static_cast<u128>(a + c) +
                                 static_cast<u128>(b + c) * static_cast<u128>(b + c);
                const u128 rhs =
                    4U * static_cast<u128>(b + c) * static_cast<u128>(a + c);
                if (lhs == rhs) {
                    total += static_cast<u128>(a + b + c);
                }
            }
        }
    }
    return total;
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string digits;
    while (value != 0) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        digits.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(digits.begin(), digits.end());
    return digits;
}

bool run_checkpoints() {
    if (solve(50ULL, 1U) != brute_force(50ULL)) {
        std::cerr << "Checkpoint failed for limit=50" << '\n';
        return false;
    }
    if (solve(200ULL, 1U) != brute_force(200ULL)) {
        std::cerr << "Checkpoint failed for limit=200" << '\n';
        return false;
    }
    if (solve(10'000'000ULL, 1U) != static_cast<u128>(23'871'972'654'940ULL)) {
        std::cerr << "Checkpoint failed for limit=10000000" << '\n';
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

    const u64 max_n = isqrt_u64(options.limit) + 2ULL;
    const unsigned thread_count = resolve_thread_count(options, max_n);

    if (options.run_checkpoints && !run_checkpoints()) {
        return 2;
    }

    std::cout << to_string_u128(solve(options.limit, thread_count)) << '\n';
    return 0;
}
