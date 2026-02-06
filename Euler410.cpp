#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u8 = std::uint8_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultR1 = 100'000'000ULL;
constexpr u64 kDefaultX1 = 1'000'000'000ULL;
constexpr u64 kDefaultR2 = 1'000'000'000ULL;
constexpr u64 kDefaultX2 = 100'000'000ULL;

constexpr u64 kCheckpointR1 = 1ULL;
constexpr u64 kCheckpointX1 = 5ULL;
constexpr u64 kCheckpointExpected1 = 10ULL;
constexpr u64 kCheckpointR2 = 2ULL;
constexpr u64 kCheckpointX2 = 10ULL;
constexpr u64 kCheckpointExpected2 = 52ULL;
constexpr u64 kCheckpointR3 = 10ULL;
constexpr u64 kCheckpointX3 = 100ULL;
constexpr u64 kCheckpointExpected3 = 3384ULL;

constexpr u64 kBruteCheckR = 12ULL;
constexpr u64 kBruteCheckX = 20ULL;
constexpr u64 kBruteCheckExpected = 824ULL;

constexpr u64 kThreadConsistencyR = 200'000ULL;
constexpr u64 kThreadConsistencyX = 300'000ULL;

struct Options {
    u64 r1 = kDefaultR1;
    u64 x1 = kDefaultX1;
    u64 r2 = kDefaultR2;
    u64 x2 = kDefaultX2;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
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

bool parse_arguments(int argc, char** argv, Options& options) {
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
        if (parse_u64_after_prefix(arg, "--r1=", parsed_u64)) {
            options.r1 = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--x1=", parsed_u64)) {
            options.x1 = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--r2=", parsed_u64)) {
            options.r2 = parsed_u64;
            continue;
        }
        if (parse_u64_after_prefix(arg, "--x2=", parsed_u64)) {
            options.x2 = parsed_u64;
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

    if (options.r1 == 0ULL || options.x1 == 0ULL || options.r2 == 0ULL || options.x2 == 0ULL) {
        std::cerr << "All bounds must be at least 1.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2ULL) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads == 0U) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0U) {
            threads = 1U;
        }
    }

    return std::max(1U, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
}

std::string to_string_u128(u128 value) {
    if (value == 0) {
        return "0";
    }

    std::string s;
    while (value > 0) {
        const unsigned digit = static_cast<unsigned>(value % 10U);
        s.push_back(static_cast<char>('0' + digit));
        value /= 10U;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

std::vector<u8> build_distinct_prime_counts(u64 limit) {
    std::vector<u8> omega(static_cast<std::size_t>(limit) + 1ULL, 0U);
    if (limit < 2ULL) {
        return omega;
    }

    std::vector<bool> composite((static_cast<std::size_t>(limit) >> 1ULL) + 1ULL, false);
    const u64 root = static_cast<u64>(std::sqrt(static_cast<long double>(limit)));

    for (u64 p = 3ULL; p <= root; p += 2ULL) {
        if (composite[static_cast<std::size_t>(p >> 1ULL)]) {
            continue;
        }

        const u64 step = p << 1ULL;
        for (u64 j = p * p; j <= limit; j += step) {
            composite[static_cast<std::size_t>(j >> 1ULL)] = true;
        }
    }

    for (u64 m = 2ULL; m <= limit; m += 2ULL) {
        ++omega[static_cast<std::size_t>(m)];
    }
    for (u64 p = 3ULL; p <= limit; p += 2ULL) {
        if (composite[static_cast<std::size_t>(p >> 1ULL)]) {
            continue;
        }
        for (u64 m = p; m <= limit; m += p) {
            ++omega[static_cast<std::size_t>(m)];
        }
    }

    return omega;
}

u64 brute_force_reference(u64 R, u64 X) {
    u64 total = 0ULL;

    for (u64 r = 1ULL; r <= R; ++r) {
        for (u64 a = 1ULL; a <= X; ++a) {
            const u64 a2 = a * a;

            for (u64 u = 1ULL; u <= a2; ++u) {
                if (a2 % u != 0ULL) {
                    continue;
                }
                const u64 v = a2 / u;
                const long long d = static_cast<long long>(v) - static_cast<long long>(u);
                const u64 k = u + v;

                const u64 rk = r * k;
                if (rk % a != 0ULL) {
                    continue;
                }

                const long long s_abs = static_cast<long long>(rk / a);
                if (((s_abs - d) & 1LL) != 0LL) {
                    continue;
                }

                total += 2ULL;  // s = +s_abs and s = -s_abs
            }
        }
    }

    return total;
}

u128 solve_F_with_omega(u64 R,
                        u64 X,
                        const std::vector<u8>& omega,
                        bool allow_multithreading,
                        unsigned requested_threads) {
    const u64 limit = std::min(R, X);
    if (limit == 0ULL) {
        return 0;
    }

    const unsigned threads = choose_thread_count(
        allow_multithreading, requested_threads, static_cast<std::size_t>(limit));
    std::vector<u128> partial(threads, 0);

    auto worker = [&](unsigned tid) {
        u128 local = 0;

        for (u64 m = 1ULL + static_cast<u64>(tid); m <= limit; m += threads) {
            const u64 t = R / m;
            const u64 g = X / m;
            const u64 weight = 1ULL << omega[static_cast<std::size_t>(m)];

            u64 unit = 0ULL;
            if ((m & 1ULL) != 0ULL) {
                unit = 2ULL * t * g;
            } else {
                unit = t * g + ((t & 1ULL) & (g & 1ULL));
            }

            local += static_cast<u128>(weight) * static_cast<u128>(unit);
        }

        partial[tid] = local;
    };

    std::vector<std::thread> pool;
    pool.reserve(threads > 0U ? threads - 1U : 0U);
    for (unsigned t = 1U; t < threads; ++t) {
        pool.emplace_back(worker, t);
    }
    worker(0U);
    for (std::thread& thread : pool) {
        thread.join();
    }

    u128 total = 0;
    for (const u128 value : partial) {
        total += value;
    }
    return total;
}

u128 solve_F(u64 R,
             u64 X,
             bool allow_multithreading,
             unsigned requested_threads) {
    const u64 limit = std::min(R, X);
    const std::vector<u8> omega = build_distinct_prime_counts(limit);
    return solve_F_with_omega(R, X, omega, allow_multithreading, requested_threads);
}

bool run_checkpoints(const Options& options) {
    {
        const u128 got = solve_F(kCheckpointR1,
                                 kCheckpointX1,
                                 options.allow_multithreading,
                                 options.requested_threads);
        if (got != static_cast<u128>(kCheckpointExpected1)) {
            std::cerr << "Checkpoint failed: F(" << kCheckpointR1 << ", " << kCheckpointX1
                      << ") expected " << kCheckpointExpected1
                      << ", got " << to_string_u128(got) << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: F(" << kCheckpointR1 << ", " << kCheckpointX1
                  << ") = " << kCheckpointExpected1 << '\n';
    }

    {
        const u128 got = solve_F(kCheckpointR2,
                                 kCheckpointX2,
                                 options.allow_multithreading,
                                 options.requested_threads);
        if (got != static_cast<u128>(kCheckpointExpected2)) {
            std::cerr << "Checkpoint failed: F(" << kCheckpointR2 << ", " << kCheckpointX2
                      << ") expected " << kCheckpointExpected2
                      << ", got " << to_string_u128(got) << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: F(" << kCheckpointR2 << ", " << kCheckpointX2
                  << ") = " << kCheckpointExpected2 << '\n';
    }

    {
        const u128 got = solve_F(kCheckpointR3,
                                 kCheckpointX3,
                                 options.allow_multithreading,
                                 options.requested_threads);
        if (got != static_cast<u128>(kCheckpointExpected3)) {
            std::cerr << "Checkpoint failed: F(" << kCheckpointR3 << ", " << kCheckpointX3
                      << ") expected " << kCheckpointExpected3
                      << ", got " << to_string_u128(got) << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: F(" << kCheckpointR3 << ", " << kCheckpointX3
                  << ") = " << kCheckpointExpected3 << '\n';
    }

    {
        const u64 brute = brute_force_reference(kBruteCheckR, kBruteCheckX);
        if (brute != kBruteCheckExpected) {
            std::cerr << "Internal brute checkpoint mismatch: expected "
                      << kBruteCheckExpected << ", got " << brute << '\n';
            return false;
        }

        const u128 fast = solve_F(kBruteCheckR,
                                  kBruteCheckX,
                                  options.allow_multithreading,
                                  options.requested_threads);
        if (fast != static_cast<u128>(brute)) {
            std::cerr << "Fast-vs-brute checkpoint failed at F(" << kBruteCheckR << ", "
                      << kBruteCheckX << "): brute=" << brute
                      << ", fast=" << to_string_u128(fast) << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: F(" << kBruteCheckR << ", " << kBruteCheckX
                  << ") = " << brute << " (matches brute force)\n";
    }

    if (options.allow_multithreading) {
        const u64 limit = std::min(kThreadConsistencyR, kThreadConsistencyX);
        const std::vector<u8> omega = build_distinct_prime_counts(limit);

        const u128 multi = solve_F_with_omega(kThreadConsistencyR,
                                              kThreadConsistencyX,
                                              omega,
                                              true,
                                              options.requested_threads);
        const u128 single = solve_F_with_omega(kThreadConsistencyR,
                                               kThreadConsistencyX,
                                               omega,
                                               false,
                                               options.requested_threads);
        if (multi != single) {
            std::cerr << "Thread consistency failed at F(" << kThreadConsistencyR << ", "
                      << kThreadConsistencyX << "): single=" << to_string_u128(single)
                      << ", multi=" << to_string_u128(multi) << '\n';
            return false;
        }
        std::cout << "Checkpoint passed: thread consistency at F(" << kThreadConsistencyR
                  << ", " << kThreadConsistencyX << ")\n";
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    try {
        const auto start_time = std::chrono::steady_clock::now();

        if (options.run_checkpoints && !run_checkpoints(options)) {
            return 1;
        }

        const u64 limit1 = std::min(options.r1, options.x1);
        const u64 limit2 = std::min(options.r2, options.x2);
        const u64 max_limit = std::max(limit1, limit2);
        const std::vector<u8> omega = build_distinct_prime_counts(max_limit);

        const u128 f1 = solve_F_with_omega(options.r1,
                                           options.x1,
                                           omega,
                                           options.allow_multithreading,
                                           options.requested_threads);

        u128 f2 = 0;
        if (options.r2 == options.x1 && options.x2 == options.r1) {
            f2 = f1;
        } else {
            f2 = solve_F_with_omega(options.r2,
                                    options.x2,
                                    omega,
                                    options.allow_multithreading,
                                    options.requested_threads);
        }

        const u128 answer = f1 + f2;

        const auto end_time = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end_time - start_time;

        std::cout << "F(" << options.r1 << ", " << options.x1 << ") = "
                  << to_string_u128(f1) << '\n';
        std::cout << "F(" << options.r2 << ", " << options.x2 << ") = "
                  << to_string_u128(f2) << '\n';
        std::cout << "Answer = " << to_string_u128(answer) << '\n';
        std::cout << "Elapsed: " << elapsed.count() << " seconds\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
