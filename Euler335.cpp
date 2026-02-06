#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 kExponentLimit = 1'000'000'000'000'000'000ULL;
constexpr u64 kModulus = 40'353'607ULL;  // 7^9

struct Options {
    bool run_checks = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct ValidationTask {
    u64 n = 0;
    u64 expected = 0;
    std::string label;
};

bool parse_u64_after_prefix(const std::string& arg, const char* prefix, u64& value) {
    const std::string p(prefix);
    if (arg.rfind(p, 0) != 0U) return false;

    const std::string tail = arg.substr(p.size());
    if (tail.empty()) return false;

    u64 parsed = 0;
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
    u64 parsed = 0;
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

u64 pow_u64(u64 base, unsigned exp) {
    u64 out = 1ULL;
    for (unsigned i = 0U; i < exp; ++i) out *= base;
    return out;
}

u64 simulate_M(const u64 bowls_count) {
    std::vector<u64> bowls(static_cast<std::size_t>(bowls_count), 1ULL);
    u64 pos = 0ULL;

    for (u64 steps = 1ULL;; ++steps) {
        const u64 beans = bowls[static_cast<std::size_t>(pos)];
        bowls[static_cast<std::size_t>(pos)] = 0ULL;

        for (u64 i = 0ULL; i < beans; ++i) {
            ++pos;
            if (pos == bowls_count) pos = 0ULL;
            ++bowls[static_cast<std::size_t>(pos)];
        }

        bool all_ones = true;
        for (const u64 value : bowls) {
            if (value != 1ULL) {
                all_ones = false;
                break;
            }
        }
        if (all_ones) return steps;
    }
}

u64 closed_form_M_for_power_of_two_plus_one(unsigned k) {
    // M(2^k + 1) = 2^(k+1) + 4^k - 3^k
    return pow_u64(2ULL, k + 1U) + pow_u64(4ULL, k) - pow_u64(3ULL, k);
}

bool run_validations(const bool allow_multithreading, const unsigned requested_threads) {
    std::vector<ValidationTask> tasks;
    tasks.reserve(12);

    tasks.push_back(ValidationTask{2ULL, 2ULL, "statement: M(2)"});
    tasks.push_back(ValidationTask{5ULL, 15ULL, "statement: M(5)"});
    tasks.push_back(ValidationTask{100ULL, 10'920ULL, "statement: M(100)"});

    for (unsigned k = 0U; k <= 8U; ++k) {
        const u64 n = (1ULL << k) + 1ULL;
        tasks.push_back(ValidationTask{
            n,
            closed_form_M_for_power_of_two_plus_one(k),
            "closed form at k=" + std::to_string(k),
        });
    }

    const unsigned thread_count =
        choose_thread_count(allow_multithreading, requested_threads, tasks.size());

    std::atomic<std::size_t> next_task{0ULL};
    std::vector<std::string> failures;
    std::mutex failure_mutex;
    std::vector<std::thread> pool;
    pool.reserve(thread_count);

    for (unsigned t = 0U; t < thread_count; ++t) {
        pool.emplace_back([&]() {
            while (true) {
                const std::size_t id = next_task.fetch_add(1ULL, std::memory_order_relaxed);
                if (id >= tasks.size()) break;

                const ValidationTask& task = tasks[id];
                const u64 got = simulate_M(task.n);
                if (got == task.expected) continue;

                std::ostringstream oss;
                oss << "Validation failed (" << task.label << "): M(" << task.n << ") = " << got
                    << ", expected " << task.expected;

                std::lock_guard<std::mutex> lock(failure_mutex);
                failures.push_back(oss.str());
            }
        });
    }

    for (std::thread& th : pool) th.join();

    if (failures.empty()) return true;
    std::sort(failures.begin(), failures.end());
    for (const std::string& line : failures) std::cerr << line << '\n';
    return false;
}

u64 mul_mod(const u64 a, const u64 b, const u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, const u64 mod) {
    u64 out = 1ULL % mod;
    base %= mod;

    while (exp > 0ULL) {
        if ((exp & 1ULL) != 0ULL) out = mul_mod(out, base, mod);
        base = mul_mod(base, base, mod);
        exp >>= 1ULL;
    }

    return out;
}

u64 mod_inverse(const u64 a, const u64 mod) {
    std::int64_t t = 0;
    std::int64_t new_t = 1;
    std::int64_t r = static_cast<std::int64_t>(mod);
    std::int64_t new_r = static_cast<std::int64_t>(a % mod);

    while (new_r != 0) {
        const std::int64_t q = r / new_r;
        const std::int64_t temp_t = t - q * new_t;
        t = new_t;
        new_t = temp_t;

        const std::int64_t temp_r = r - q * new_r;
        r = new_r;
        new_r = temp_r;
    }

    if (r != 1) return 0;
    if (t < 0) t += static_cast<std::int64_t>(mod);
    return static_cast<u64>(t);
}

u64 solve() {
    const u64 inv2 = mod_inverse(2ULL, kModulus);
    const u64 inv3 = mod_inverse(3ULL, kModulus);

    const u64 p2 = pow_mod(2ULL, kExponentLimit + 1ULL, kModulus);
    const u64 p3 = pow_mod(3ULL, kExponentLimit + 1ULL, kModulus);
    const u64 p4 = pow_mod(4ULL, kExponentLimit + 1ULL, kModulus);

    const u64 sum_2 = mul_mod(2ULL, (p2 + kModulus - 1ULL) % kModulus, kModulus);
    const u64 sum_4 = mul_mod((p4 + kModulus - 1ULL) % kModulus, inv3, kModulus);
    const u64 sum_3 = mul_mod((p3 + kModulus - 1ULL) % kModulus, inv2, kModulus);

    return (sum_2 + sum_4 + kModulus - sum_3) % kModulus;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) return 1;

    if (options.run_checks && !run_validations(options.allow_multithreading, options.requested_threads)) {
        return 1;
    }

    std::cout << solve() << '\n';
    return 0;
}
