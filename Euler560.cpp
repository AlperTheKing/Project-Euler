#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 kMod = 1'000'000'007U;
constexpr u32 kDefaultN = 10'000'000U;
constexpr u32 kDefaultK = 10'000'000U;

constexpr u32 kCheckpointN1 = 5U;
constexpr u32 kCheckpointK1 = 2U;
constexpr u32 kCheckpointExpected1 = 6U;

constexpr u32 kCheckpointN2 = 10U;
constexpr u32 kCheckpointK2 = 5U;
constexpr u32 kCheckpointExpected2 = 9'964U;

constexpr u32 kCheckpointN3 = 10U;
constexpr u32 kCheckpointK3 = 10U;
constexpr u32 kCheckpointExpected3 = 472'400'303U;

constexpr u32 kCheckpointN4 = 1'000U;
constexpr u32 kCheckpointK4 = 1'000U;
constexpr u32 kCheckpointExpected4 = 954'021'836U;

constexpr u32 kFormulaValidationLimit = 180U;
constexpr u32 kBruteCrossN = 30U;
constexpr u32 kBruteCrossK = 4U;
constexpr u32 kThreadConsistencyN = 100'000U;
constexpr u32 kThreadConsistencyK = 100'000U;

struct Options {
    u32 n = kDefaultN;
    u32 k = kDefaultK;
    bool run_checkpoints = true;
    bool allow_multithreading = true;
    unsigned requested_threads = 0U;
};

struct SieveData {
    u32 limit = 0U;
    std::vector<u32> spf;
    std::vector<u32> pi;
};

inline u32 add_mod(const u32 a, const u32 b) {
    const u32 s = a + b;
    return (s >= kMod) ? (s - kMod) : s;
}

inline u32 sub_mod(const u32 a, const u32 b) {
    return (a >= b) ? (a - b) : (a + kMod - b);
}

inline u32 mul_mod(const u32 a, const u32 b) {
    return static_cast<u32>((static_cast<u64>(a) * static_cast<u64>(b)) % kMod);
}

u32 mod_pow(u32 base, u32 exponent) {
    u32 result = 1U;
    u32 b = base;
    u32 e = exponent;

    while (e > 0U) {
        if ((e & 1U) != 0U) {
            result = mul_mod(result, b);
        }
        b = mul_mod(b, b);
        e >>= 1U;
    }

    return result;
}

bool parse_u32_after_prefix(const std::string& arg, const char* prefix, u32& value) {
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
        const u32 digit = static_cast<u32>(c - '0');
        parsed = parsed * 10ULL + static_cast<u64>(digit);
        if (parsed > static_cast<u64>(std::numeric_limits<u32>::max())) {
            return false;
        }
    }

    value = static_cast<u32>(parsed);
    return true;
}

bool parse_unsigned_after_prefix(const std::string& arg,
                                 const char* prefix,
                                 unsigned& value) {
    u32 parsed = 0U;
    if (!parse_u32_after_prefix(arg, prefix, parsed)) {
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

        u32 parsed_u32 = 0U;
        if (parse_u32_after_prefix(arg, "--n=", parsed_u32)) {
            options.n = parsed_u32;
            continue;
        }
        if (parse_u32_after_prefix(arg, "--k=", parsed_u32)) {
            options.k = parsed_u32;
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

    if (options.n < 2U) {
        std::cerr << "--n must be at least 2.\n";
        return false;
    }
    if (options.k == 0U) {
        std::cerr << "--k must be at least 1.\n";
        return false;
    }

    return true;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 8'192ULL) {
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

SieveData build_linear_sieve(const u32 limit) {
    SieveData data;
    data.limit = limit;
    data.spf.assign(static_cast<std::size_t>(limit) + 1ULL, 0U);
    data.pi.assign(static_cast<std::size_t>(limit) + 1ULL, 0U);

    if (limit < 2U) {
        return data;
    }

    std::vector<u32> primes;
    primes.reserve(static_cast<std::size_t>(limit / 10U) + 32ULL);

    for (u32 i = 2U; i <= limit; ++i) {
        if (data.spf[i] == 0U) {
            data.spf[i] = i;
            primes.push_back(i);
        }

        for (const u32 p : primes) {
            const u64 v = static_cast<u64>(i) * static_cast<u64>(p);
            if (v > static_cast<u64>(limit)) {
                break;
            }
            data.spf[static_cast<std::size_t>(v)] = p;
            if (p == data.spf[i]) {
                break;
            }
        }
    }

    u32 prime_count = 0U;
    for (u32 x = 1U; x <= limit; ++x) {
        if (x >= 2U && data.spf[x] == x) {
            ++prime_count;
        }
        data.pi[x] = prime_count;
    }

    return data;
}

u32 grundy_formula(const u32 x, const SieveData& sieve) {
    if (x == 1U) {
        return 1U;
    }
    if ((x & 1U) == 0U) {
        return 0U;
    }
    return sieve.pi[sieve.spf[x]];
}

std::vector<u32> build_grundy_frequencies(const u32 n) {
    const u32 limit = n - 1U;
    const SieveData sieve = build_linear_sieve(limit);

    const u32 max_index = std::max(1U, sieve.pi[limit]);
    std::vector<u32> counts(static_cast<std::size_t>(max_index) + 1ULL, 0U);

    counts[0] = limit / 2U;
    counts[1] = 1U;

    for (u32 x = 3U; x <= limit; x += 2U) {
        const u32 p = sieve.spf[x];
        const u32 idx = sieve.pi[p];
        ++counts[idx];
    }

    return counts;
}

void fwht_xor(std::vector<u32>& a) {
    const std::size_t n = a.size();

    for (std::size_t len = 1ULL; len < n; len <<= 1ULL) {
        const std::size_t step = len << 1ULL;
        for (std::size_t i = 0ULL; i < n; i += step) {
            for (std::size_t j = 0ULL; j < len; ++j) {
                const u32 u = a[i + j];
                const u32 v = a[i + j + len];
                a[i + j] = add_mod(u, v);
                a[i + j + len] = sub_mod(u, v);
            }
        }
    }
}

void pow_coefficients(std::vector<u32>& coeffs,
                      const u32 exponent,
                      const bool allow_multithreading,
                      const unsigned requested_threads) {
    const std::size_t n = coeffs.size();
    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, n);

    if (threads == 1U) {
        for (u32& v : coeffs) {
            v = mod_pow(v, exponent);
        }
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(threads);

    for (unsigned t = 0U; t < threads; ++t) {
        const std::size_t begin = (n * static_cast<std::size_t>(t)) / static_cast<std::size_t>(threads);
        const std::size_t end = (n * static_cast<std::size_t>(t + 1U)) /
                                static_cast<std::size_t>(threads);

        workers.emplace_back([&coeffs, exponent, begin, end]() {
            for (std::size_t i = begin; i < end; ++i) {
                coeffs[i] = mod_pow(coeffs[i], exponent);
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }
}

u32 solve_losing_positions(const u32 n,
                          const u32 k,
                          const bool allow_multithreading,
                          const unsigned requested_threads) {
    const std::vector<u32> counts = build_grundy_frequencies(n);

    std::size_t size = 1ULL;
    while (size < counts.size()) {
        size <<= 1ULL;
    }

    std::vector<u32> transformed(size, 0U);
    for (std::size_t i = 0ULL; i < counts.size(); ++i) {
        transformed[i] = counts[i] % kMod;
    }

    fwht_xor(transformed);
    pow_coefficients(transformed, k, allow_multithreading, requested_threads);
    fwht_xor(transformed);

    const u32 inv_size = mod_pow(static_cast<u32>(size % static_cast<std::size_t>(kMod)), kMod - 2U);
    return mul_mod(transformed[0], inv_size);
}

std::vector<u32> brute_grundy_values(const u32 limit) {
    std::vector<u32> grundy(static_cast<std::size_t>(limit) + 1ULL, 0U);

    for (u32 n = 1U; n <= limit; ++n) {
        std::vector<unsigned char> seen(static_cast<std::size_t>(n) + 8ULL, 0U);

        for (u32 take = 1U; take <= n; ++take) {
            if (std::gcd(take, n) != 1U) {
                continue;
            }
            const u32 next_value = n - take;
            const u32 g = grundy[next_value];
            if (g >= seen.size()) {
                seen.resize(static_cast<std::size_t>(g) + 1ULL, 0U);
            }
            seen[g] = 1U;
        }

        u32 mex = 0U;
        while (mex < seen.size() && seen[mex] != 0U) {
            ++mex;
        }
        grundy[n] = mex;
    }

    return grundy;
}

u32 brute_losing_positions_small(const u32 n, const u32 k) {
    const u32 limit = n - 1U;
    const std::vector<u32> grundy = brute_grundy_values(limit);

    u32 max_grundy = 0U;
    for (u32 x = 1U; x <= limit; ++x) {
        max_grundy = std::max(max_grundy, grundy[x]);
    }

    std::size_t size = 1ULL;
    while (size <= static_cast<std::size_t>(max_grundy)) {
        size <<= 1ULL;
    }

    std::vector<u32> counts(size, 0U);
    for (u32 x = 1U; x <= limit; ++x) {
        ++counts[grundy[x]];
    }

    std::vector<u32> dp(size, 0U);
    std::vector<u32> next(size, 0U);
    dp[0] = 1U;

    for (u32 pile = 0U; pile < k; ++pile) {
        std::fill(next.begin(), next.end(), 0U);
        for (std::size_t xor_value = 0ULL; xor_value < size; ++xor_value) {
            const u32 ways = dp[xor_value];
            if (ways == 0U) {
                continue;
            }
            for (std::size_t g = 0ULL; g < size; ++g) {
                const u32 freq = counts[g];
                if (freq == 0U) {
                    continue;
                }
                const std::size_t target = xor_value ^ g;
                const u32 contribution = mul_mod(ways, freq % kMod);
                next[target] = add_mod(next[target], contribution);
            }
        }
        dp.swap(next);
    }

    return dp[0];
}

bool expect_equal(const char* label, const u32 actual, const u32 expected) {
    if (actual == expected) {
        return true;
    }

    std::cerr << "Checkpoint failed for " << label << ": got " << actual
              << ", expected " << expected << '\n';
    return false;
}

bool run_checkpoints() {
    bool ok = true;

    {
        const SieveData sieve = build_linear_sieve(kFormulaValidationLimit);
        const std::vector<u32> brute = brute_grundy_values(kFormulaValidationLimit);

        for (u32 x = 1U; x <= kFormulaValidationLimit; ++x) {
            const u32 predicted = grundy_formula(x, sieve);
            if (predicted != brute[x]) {
                std::cerr << "Grundy formula mismatch at n=" << x << ": got " << predicted
                          << ", expected " << brute[x] << '\n';
                ok = false;
                break;
            }
        }
    }

    ok = ok && expect_equal("L(5,2)",
                            solve_losing_positions(kCheckpointN1, kCheckpointK1, false, 1U),
                            kCheckpointExpected1);
    ok = ok && expect_equal("L(10,5)",
                            solve_losing_positions(kCheckpointN2, kCheckpointK2, false, 1U),
                            kCheckpointExpected2);
    ok = ok && expect_equal("L(10,10)",
                            solve_losing_positions(kCheckpointN3, kCheckpointK3, false, 1U),
                            kCheckpointExpected3);
    ok = ok && expect_equal("L(10^3,10^3)",
                            solve_losing_positions(kCheckpointN4, kCheckpointK4, false, 1U),
                            kCheckpointExpected4);

    const u32 brute_cross = brute_losing_positions_small(kBruteCrossN, kBruteCrossK);
    const u32 fast_cross = solve_losing_positions(kBruteCrossN, kBruteCrossK, false, 1U);
    ok = ok && expect_equal("brute-vs-fast", fast_cross, brute_cross);

    unsigned hw_threads = std::thread::hardware_concurrency();
    if (hw_threads == 0U) {
        hw_threads = 1U;
    }
    if (hw_threads >= 2U) {
        const u32 single =
            solve_losing_positions(kThreadConsistencyN, kThreadConsistencyK, false, 1U);
        const u32 threaded =
            solve_losing_positions(kThreadConsistencyN, kThreadConsistencyK, true, 2U);
        ok = ok && expect_equal("thread-consistency", threaded, single);
    }

    if (!ok) {
        std::cerr << "At least one checkpoint failed.\n";
    }
    return ok;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }

    const auto start = std::chrono::steady_clock::now();
    const u32 answer = solve_losing_positions(
        options.n, options.k, options.allow_multithreading, options.requested_threads);
    const auto end = std::chrono::steady_clock::now();

    const std::chrono::duration<long double> elapsed = end - start;
    std::cout << "L(" << options.n << ", " << options.k << ") mod " << kMod << " = "
              << answer << '\n';
    std::cout << "Elapsed: " << elapsed.count() << " seconds\n";

    return 0;
}
