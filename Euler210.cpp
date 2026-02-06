#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;

constexpr u64 kDefaultR = 1'000'000'000ULL;
constexpr u64 kCheckpointR1 = 4ULL;
constexpr u64 kCheckpointExpected1 = 24ULL;
constexpr u64 kCheckpointR2 = 8ULL;
constexpr u64 kCheckpointExpected2 = 100ULL;
constexpr u64 kThreadConsistencyR = 4'000'000ULL;

struct Options {
    u64 r = kDefaultR;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0) {
        return false;
    }

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) {
        return false;
    }

    u64 parsed = 0;
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
    u64 parsed = 0;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }

    if (parsed > static_cast<u64>(std::numeric_limits<unsigned>::max())) {
        return false;
    }

    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_arguments(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);

        if (arg == "--single-thread") {
            options.allow_multithreading = false;
            continue;
        }
        if (arg == "--skip-checkpoints") {
            options.run_checkpoints = false;
            continue;
        }

        u64 r = 0;
        if (parse_u64_after_prefix(arg, "--r=", r)) {
            options.r = r;
            continue;
        }

        unsigned threads = 0;
        if (parse_unsigned_after_prefix(arg, "--threads=", threads)) {
            options.requested_threads = threads;
            continue;
        }

        std::cerr << "Unknown argument: " << arg << '\n';
        return false;
    }

    if (options.r == 0ULL) {
        std::cerr << "--r must be >= 1.\n";
        return false;
    }
    if ((options.r % 4ULL) != 0ULL) {
        std::cerr << "This solver requires r to be divisible by 4.\n";
        return false;
    }

    return true;
}

u64 isqrt_u64(u64 n) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while ((r + 1ULL) <= n / (r + 1ULL)) {
        ++r;
    }
    while (r > n / r) {
        --r;
    }
    return r;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2'000'000ULL) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    return std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

i128 count_circle_u_range(u64 limit, int parity, u64 begin_u, u64 end_u) {
    if (begin_u > end_u) {
        return 0;
    }

    constexpr u64 kNoValue = std::numeric_limits<u64>::max();

    u64 v = isqrt_u64(limit - begin_u * begin_u);
    if (static_cast<int>(v & 1ULL) != parity) {
        if (v == 0ULL) {
            v = kNoValue;
        } else {
            --v;
        }
    }

    i128 total = 0;
    for (u64 u = begin_u; u <= end_u; u += 2ULL) {
        if (v == kNoValue) {
            break;
        }

        const u64 u2 = u * u;
        while (u2 + v * v > limit) {
            if (v <= 1ULL) {
                v = kNoValue;
                break;
            }
            v -= 2ULL;
        }
        if (v == kNoValue) {
            break;
        }

        const i128 count_v = static_cast<i128>(v + 1ULL);
        const i128 multiplicity = (u == 0ULL ? static_cast<i128>(1) : static_cast<i128>(2));
        total += multiplicity * count_v;
    }

    return total;
}

u64 count_circle_points_with_parity(u64 a,
                                    bool allow_multithreading,
                                    unsigned requested_threads) {
    const u64 limit = 2ULL * a * a - 1ULL;
    const int parity = static_cast<int>(a & 1ULL);

    u64 max_u = isqrt_u64(limit);
    if (static_cast<int>(max_u & 1ULL) != parity) {
        --max_u;
    }

    const u64 first_u = static_cast<u64>(parity);
    if (first_u > max_u) {
        return 0ULL;
    }

    const u64 step_count = (max_u - first_u) / 2ULL + 1ULL;
    const unsigned threads =
        choose_thread_count(allow_multithreading, requested_threads, static_cast<std::size_t>(step_count));

    i128 total = 0;
    if (threads == 1U) {
        total = count_circle_u_range(limit, parity, first_u, max_u);
    } else {
        std::vector<std::thread> workers;
        std::vector<i128> partial(threads, 0);
        workers.reserve(threads);

        for (unsigned t = 0; t < threads; ++t) {
            const u64 step_begin = static_cast<u64>((static_cast<__int128>(step_count) * t) / threads);
            const u64 step_end_exclusive =
                static_cast<u64>((static_cast<__int128>(step_count) * (t + 1ULL)) / threads);

            if (step_begin >= step_end_exclusive) {
                continue;
            }

            const u64 begin_u = first_u + 2ULL * step_begin;
            const u64 end_u = first_u + 2ULL * (step_end_exclusive - 1ULL);
            workers.emplace_back([&, t, begin_u, end_u]() {
                partial[t] = count_circle_u_range(limit, parity, begin_u, end_u);
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        for (const i128 value : partial) {
            total += value;
        }
    }

    return static_cast<u64>(total);
}

i128 solve_obtuse_count(u64 r,
                        bool allow_multithreading,
                        unsigned requested_threads) {
    const u64 a = r / 4ULL;

    // Angle at O is obtuse when x + y < 0.
    const i128 count_o = static_cast<i128>(r) * static_cast<i128>(r) +
                         static_cast<i128>(r / 2ULL);

    // Angle at C is obtuse when x + y > r/2.
    const i128 count_c =
        static_cast<i128>(r) * static_cast<i128>(2ULL * r + 1ULL) / static_cast<i128>(4);

    // Angle at B is obtuse when B is strictly inside the circle with diameter OC.
    const i128 count_b = static_cast<i128>(count_circle_points_with_parity(
        a,
        allow_multithreading,
        requested_threads));

    // Degenerate collinear points on y = x were included once above but must be excluded.
    const i128 degenerate = static_cast<i128>(r - 1ULL);

    return count_o + count_c + count_b - degenerate;
}

i128 brute_force_count(u64 r) {
    const i64 a = static_cast<i64>(r / 4ULL);

    i128 count = 0;
    for (i64 x = -static_cast<i64>(r); x <= static_cast<i64>(r); ++x) {
        for (i64 y = -static_cast<i64>(r); y <= static_cast<i64>(r); ++y) {
            if (std::llabs(x) + std::llabs(y) > static_cast<i64>(r)) {
                continue;
            }

            if (x == y) {
                // O, B, C are collinear -> largest angle is 180 degrees, excluded.
                continue;
            }

            const i64 sum = x + y;
            const i64 dot_o = a * sum;
            const i64 dot_b = x * x + y * y - a * sum;
            const i64 dot_c = 2LL * a * a - a * sum;

            if (dot_o < 0LL || dot_b < 0LL || dot_c < 0LL) {
                ++count;
            }
        }
    }

    return count;
}

std::string to_string_i128(i128 value) {
    if (value == 0) {
        return "0";
    }

    bool negative = value < 0;
    if (negative) {
        value = -value;
    }

    std::string out;
    while (value > 0) {
        const int digit = static_cast<int>(value % 10);
        out.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }

    if (negative) {
        out.push_back('-');
    }
    std::reverse(out.begin(), out.end());
    return out;
}

bool run_checkpoints(const Options& options) {
    struct FixedCheckpoint {
        u64 r;
        u64 expected;
    };

    const std::vector<FixedCheckpoint> fixed = {
        {kCheckpointR1, kCheckpointExpected1},
        {kCheckpointR2, kCheckpointExpected2},
    };

    for (const FixedCheckpoint& cp : fixed) {
        const i128 got = solve_obtuse_count(cp.r, false, 1U);
        if (got != static_cast<i128>(cp.expected)) {
            std::cerr << "Checkpoint failed: N(" << cp.r << ") expected " << cp.expected
                      << ", got " << to_string_i128(got) << '\n';
            return false;
        }
        std::cout << "Checkpoint OK: N(" << cp.r << ") = " << cp.expected << '\n';
    }

    const std::vector<u64> brute_cases = {12ULL, 16ULL, 20ULL, 24ULL, 28ULL, 32ULL};
    for (const u64 r : brute_cases) {
        const i128 brute = brute_force_count(r);
        const i128 fast = solve_obtuse_count(r, false, 1U);
        if (brute != fast) {
            std::cerr << "Brute checkpoint failed: N(" << r << ") brute=" << to_string_i128(brute)
                      << ", fast=" << to_string_i128(fast) << '\n';
            return false;
        }
        std::cout << "Checkpoint OK: brute cross-check N(" << r << ") = "
                  << to_string_i128(fast) << '\n';
    }

    if (options.allow_multithreading) {
        const i128 single_thread = solve_obtuse_count(kThreadConsistencyR, false, 1U);
        const i128 multi_thread =
            solve_obtuse_count(kThreadConsistencyR, true, options.requested_threads);
        if (single_thread != multi_thread) {
            std::cerr << "Thread-consistency checkpoint failed at N(" << kThreadConsistencyR
                      << "): single=" << to_string_i128(single_thread)
                      << ", multi=" << to_string_i128(multi_thread) << '\n';
            return false;
        }
        std::cout << "Checkpoint OK: threaded consistency at N(" << kThreadConsistencyR
                  << ") = " << to_string_i128(single_thread) << '\n';
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    const auto start = std::chrono::steady_clock::now();

    if (options.run_checkpoints) {
        if (!run_checkpoints(options)) {
            return 1;
        }
    }

    const i128 answer =
        solve_obtuse_count(options.r, options.allow_multithreading, options.requested_threads);

    const auto finish = std::chrono::steady_clock::now();
    const std::chrono::duration<long double> elapsed = finish - start;

    std::cout << "Answer: " << to_string_i128(answer) << '\n';
    std::cout << "N(" << options.r << ") = " << to_string_i128(answer) << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " s\n";

    return 0;
}
