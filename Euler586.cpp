#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <thread>
#include <vector>

namespace {

using u64 = std::uint64_t;
using u32 = std::uint32_t;
using u128 = unsigned __int128;

constexpr u64 kDefaultLimit = 1'000'000'000'000'000ULL;
constexpr int kDefaultR = 40;

constexpr u64 kCheckpointN1 = 100'000ULL;
constexpr int kCheckpointR1 = 4;
constexpr u64 kCheckpointExpected1 = 237ULL;

constexpr u64 kCheckpointN2 = 100'000'000ULL;
constexpr int kCheckpointR2 = 6;
constexpr u64 kCheckpointExpected2 = 59'517ULL;

struct Options {
    u64 limit = kDefaultLimit;
    int r = kDefaultR;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
};

struct SolveSetup {
    std::vector<std::vector<int>> exponent_patterns;
    u64 min_split_product = 0;
    int split_prime_limit = 0;
};

struct Task {
    std::vector<int> permutation;
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

bool parse_int_after_prefix(const std::string& arg, const char* prefix, int& value) {
    u64 parsed = 0;
    if (!parse_u64_after_prefix(arg, prefix, parsed)) {
        return false;
    }
    if (parsed > static_cast<u64>(std::numeric_limits<int>::max())) {
        return false;
    }
    value = static_cast<int>(parsed);
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

        u64 limit = 0;
        if (parse_u64_after_prefix(arg, "--limit=", limit)) {
            options.limit = limit;
            continue;
        }

        int r = 0;
        if (parse_int_after_prefix(arg, "--r=", r)) {
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

    if (options.r < 1) {
        std::cerr << "--r must be >= 1.\n";
        return false;
    }
    if (options.limit == 0ULL) {
        std::cerr << "--limit must be >= 1.\n";
        return false;
    }

    return true;
}

u64 isqrt_u64(u64 x) {
    u64 r = static_cast<u64>(std::sqrt(static_cast<long double>(x)));
    while ((r + 1ULL) <= x / (r + 1ULL)) {
        ++r;
    }
    while (r > x / r) {
        --r;
    }
    return r;
}

u64 pow_with_limit(u64 base, int exp, u64 limit) {
    u64 result = 1ULL;
    for (int i = 0; i < exp; ++i) {
        if (result > limit / base) {
            return limit + 1ULL;
        }
        result *= base;
    }
    return result;
}

u64 floor_nth_root_u64(u64 n, int exponent) {
    if (exponent <= 1 || n <= 1ULL) {
        return n;
    }

    u64 x = static_cast<u64>(
        std::pow(static_cast<long double>(n), 1.0L / static_cast<long double>(exponent)));
    if (x == 0ULL) {
        x = 1ULL;
    }

    while (pow_with_limit(x + 1ULL, exponent, n) <= n) {
        ++x;
    }
    while (pow_with_limit(x, exponent, n) > n) {
        --x;
    }
    return x;
}

std::vector<int> sieve_primes(int limit) {
    if (limit < 2) {
        return {};
    }

    std::vector<std::uint8_t> is_prime(static_cast<std::size_t>(limit) + 1ULL, 1U);
    is_prime[0] = 0U;
    is_prime[1] = 0U;

    const int root = static_cast<int>(std::sqrt(static_cast<long double>(limit)));
    for (int p = 2; p <= root; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] == 0U) {
            continue;
        }
        for (int64_t m = static_cast<int64_t>(p) * static_cast<int64_t>(p); m <= limit; m += p) {
            is_prime[static_cast<std::size_t>(m)] = 0U;
        }
    }

    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(limit / std::max(1.0L, std::log(static_cast<long double>(limit)))));
    for (int p = 2; p <= limit; ++p) {
        if (is_prime[static_cast<std::size_t>(p)] != 0U) {
            primes.push_back(p);
        }
    }
    return primes;
}

std::vector<int> split_primes_up_to(int limit) {
    const std::vector<int> primes = sieve_primes(limit);
    std::vector<int> split;
    split.reserve(primes.size() / 2 + 1);
    for (const int p : primes) {
        const int r = p % 5;
        if (r == 1 || r == 4) {
            split.push_back(p);
        }
    }
    return split;
}

void generate_exponent_patterns_rec(int remaining_tau,
                                    int min_factor,
                                    std::vector<int>& factors,
                                    std::set<std::vector<int>>& out_patterns) {
    if (remaining_tau == 1) {
        if (!factors.empty()) {
            std::vector<int> exponents;
            exponents.reserve(factors.size());
            for (const int f : factors) {
                exponents.push_back(f - 1);
            }
            std::sort(exponents.begin(), exponents.end(), std::greater<int>());
            out_patterns.insert(exponents);
        }
        return;
    }

    for (int f = min_factor; f <= remaining_tau; ++f) {
        if (f < 2 || (remaining_tau % f) != 0) {
            continue;
        }
        factors.push_back(f);
        generate_exponent_patterns_rec(remaining_tau / f, f, factors, out_patterns);
        factors.pop_back();
    }
}

std::vector<std::vector<int>> generate_exponent_patterns(int tau) {
    std::set<std::vector<int>> patterns;
    std::vector<int> factors;
    generate_exponent_patterns_rec(tau, 2, factors, patterns);
    return std::vector<std::vector<int>>(patterns.begin(), patterns.end());
}

bool multiply_with_limit(u64& value, u64 factor, u64 limit) {
    if (factor == 0ULL) {
        value = limit + 1ULL;
        return false;
    }
    if (value > limit / factor) {
        value = limit + 1ULL;
        return false;
    }
    value *= factor;
    return true;
}

u64 minimal_product_for_pattern(const std::vector<int>& exponents,
                                const std::vector<int>& split_seed,
                                u64 limit) {
    if (exponents.size() > split_seed.size()) {
        return limit + 1ULL;
    }

    std::vector<int> sorted_exponents = exponents;
    std::sort(sorted_exponents.begin(), sorted_exponents.end(), std::greater<int>());

    u64 product = 1ULL;
    for (std::size_t i = 0; i < sorted_exponents.size(); ++i) {
        const u64 powv = pow_with_limit(static_cast<u64>(split_seed[i]), sorted_exponents[i], limit);
        if (!multiply_with_limit(product, powv, limit)) {
            return limit + 1ULL;
        }
    }

    return product;
}

int estimate_split_prime_limit(const std::vector<std::vector<int>>& patterns,
                               u64 limit,
                               const std::vector<int>& split_seed) {
    int max_bound = 0;
    std::size_t max_needed_count = 0;

    for (const auto& pattern : patterns) {
        max_needed_count = std::max(max_needed_count, pattern.size());

        for (std::size_t j = 0; j < pattern.size(); ++j) {
            std::vector<int> others;
            others.reserve(pattern.size() - 1);
            for (std::size_t t = 0; t < pattern.size(); ++t) {
                if (t != j) {
                    others.push_back(pattern[t]);
                }
            }
            std::sort(others.begin(), others.end(), std::greater<int>());

            u64 fixed_product = 1ULL;
            bool ok = true;
            for (std::size_t t = 0; t < others.size(); ++t) {
                if (t >= split_seed.size()) {
                    ok = false;
                    break;
                }
                const u64 powv = pow_with_limit(static_cast<u64>(split_seed[t]), others[t], limit);
                if (!multiply_with_limit(fixed_product, powv, limit)) {
                    ok = false;
                    break;
                }
            }
            if (!ok || fixed_product > limit) {
                continue;
            }

            const u64 remaining = limit / fixed_product;
            const int exponent = pattern[j];
            const u64 bound = floor_nth_root_u64(remaining, exponent);
            if (bound > static_cast<u64>(std::numeric_limits<int>::max())) {
                max_bound = std::numeric_limits<int>::max();
            } else {
                max_bound = std::max(max_bound, static_cast<int>(bound));
            }
        }
    }

    if (max_needed_count > 0 && max_needed_count <= split_seed.size()) {
        max_bound = std::max(max_bound, split_seed[max_needed_count - 1]);
    }
    return max_bound;
}

SolveSetup build_setup(u64 limit, int r) {
    SolveSetup setup;

    const std::vector<int> split_seed = split_primes_up_to(500);
    if (split_seed.size() < 8) {
        return setup;
    }

    std::set<std::vector<int>> patterns;
    const std::vector<int> target_tau = {2 * r, 2 * r + 1};

    u64 min_product = std::numeric_limits<u64>::max();
    for (const int tau : target_tau) {
        const std::vector<std::vector<int>> tau_patterns = generate_exponent_patterns(tau);
        for (const auto& pattern : tau_patterns) {
            const u64 m = minimal_product_for_pattern(pattern, split_seed, limit);
            if (m <= limit) {
                patterns.insert(pattern);
                min_product = std::min(min_product, m);
            }
        }
    }

    if (patterns.empty()) {
        return setup;
    }

    setup.exponent_patterns.assign(patterns.begin(), patterns.end());
    setup.min_split_product = min_product;
    setup.split_prime_limit = estimate_split_prime_limit(setup.exponent_patterns, limit, split_seed);
    return setup;
}

std::vector<u32> build_inert_only_prefix(u64 max_y) {
    if (max_y > static_cast<u64>(std::numeric_limits<int>::max())) {
        return {};
    }
    const int y_limit = static_cast<int>(max_y);
    const std::vector<int> primes = sieve_primes(y_limit);

    std::vector<std::uint8_t> inert_only(static_cast<std::size_t>(y_limit) + 1ULL, 1U);
    inert_only[0] = 0U;

    for (const int p : primes) {
        const int mod = p % 5;
        if (p == 5 || mod == 1 || mod == 4) {
            for (int64_t m = p; m <= y_limit; m += p) {
                inert_only[static_cast<std::size_t>(m)] = 0U;
            }
        }
    }

    std::vector<u32> prefix(static_cast<std::size_t>(y_limit) + 1ULL, 0U);
    u32 running = 0U;
    for (int i = 0; i <= y_limit; ++i) {
        running += static_cast<u32>(inert_only[static_cast<std::size_t>(i)]);
        prefix[static_cast<std::size_t>(i)] = running;
    }

    return prefix;
}

std::vector<u64> build_multiplier_table(u64 max_x, const std::vector<u32>& inert_prefix) {
    std::vector<u64> g(static_cast<std::size_t>(max_x) + 1ULL, 0ULL);
    for (u64 x = 1; x <= max_x; ++x) {
        u64 total = 0ULL;
        u64 pow5 = 1ULL;
        while (pow5 <= x) {
            const u64 y = isqrt_u64(x / pow5);
            total += inert_prefix[static_cast<std::size_t>(y)];
            if (pow5 > x / 5ULL) {
                break;
            }
            pow5 *= 5ULL;
        }
        g[static_cast<std::size_t>(x)] = total;
    }
    return g;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload_hint) {
    if (!allow_multithreading) {
        return 1U;
    }

    unsigned threads = requested_threads;
    if (threads > 0U) {
        return std::max(1U, threads);
    }

    if (workload_hint < 200'000ULL) {
        return 1U;
    }

    threads = std::thread::hardware_concurrency();
    if (threads == 0U) {
        threads = 1U;
    }
    return std::max(1U, threads);
}

struct Evaluator {
    u64 n_limit = 0ULL;
    std::vector<int> split_primes;
    std::vector<int> exponent_values;
    std::vector<std::vector<u64>> power_table;
    std::vector<int> exponent_to_slot;
    std::vector<u64> multiplier_table;

    const std::vector<u64>& powers_for_exponent(int exponent) const {
        return power_table[static_cast<std::size_t>(exponent_to_slot[static_cast<std::size_t>(exponent)])];
    }

    bool tail_fits(const std::vector<int>& perm,
                   int next_pos,
                   std::size_t start_idx,
                   u64 limit) const {
        u64 product = 1ULL;
        std::size_t idx = start_idx;
        for (int pos = next_pos; pos < static_cast<int>(perm.size()); ++pos) {
            if (idx >= split_primes.size()) {
                return false;
            }
            const auto& powv = powers_for_exponent(perm[pos]);
            const u64 pe = powv[idx];
            if (pe == 0ULL || pe > limit / product) {
                return false;
            }
            product *= pe;
            ++idx;
        }
        return true;
    }

    u64 dfs_permutation(const std::vector<int>& perm,
                        int pos,
                        std::size_t start_idx,
                        u64 product) const {
        if (pos == static_cast<int>(perm.size())) {
            const u64 q = n_limit / product;
            return multiplier_table[static_cast<std::size_t>(q)];
        }

        const std::size_t remaining_slots =
            static_cast<std::size_t>(static_cast<int>(perm.size()) - pos);
        const auto& powv = powers_for_exponent(perm[pos]);
        const u64 head_limit = n_limit / product;

        u64 sum = 0ULL;
        for (std::size_t i = start_idx; i + remaining_slots <= split_primes.size(); ++i) {
            const u64 pe = powv[i];
            if (pe == 0ULL || pe > head_limit) {
                break;
            }

            const u64 new_product = product * pe;
            if (pos + 1 < static_cast<int>(perm.size())) {
                if (!tail_fits(perm, pos + 1, i + 1ULL, n_limit / new_product)) {
                    break;
                }
            }
            sum += dfs_permutation(perm, pos + 1, i + 1ULL, new_product);
        }

        return sum;
    }

    u64 evaluate_task(const Task& task) const {
        return dfs_permutation(task.permutation, 0, 0, 1ULL);
    }
};

u64 solve_f(u64 limit, int r, bool allow_multithreading, unsigned requested_threads) {
    const SolveSetup setup = build_setup(limit, r);
    if (setup.exponent_patterns.empty()) {
        return 0ULL;
    }
    if (setup.min_split_product == 0ULL || setup.min_split_product > limit) {
        return 0ULL;
    }

    const int split_prime_limit = std::max(11, setup.split_prime_limit);
    const std::vector<int> split_primes = split_primes_up_to(split_prime_limit);
    if (split_primes.empty()) {
        return 0ULL;
    }

    const u64 max_multiplier = limit / setup.min_split_product;
    const u64 max_y = isqrt_u64(max_multiplier);
    const std::vector<u32> inert_prefix = build_inert_only_prefix(max_y);
    if (inert_prefix.empty() && max_y > 0ULL) {
        return 0ULL;
    }
    const std::vector<u64> multiplier_table = build_multiplier_table(max_multiplier, inert_prefix);

    std::set<int> exponent_set;
    for (const auto& pattern : setup.exponent_patterns) {
        for (const int e : pattern) {
            exponent_set.insert(e);
        }
    }
    const std::vector<int> exponent_values(exponent_set.begin(), exponent_set.end());

    std::vector<int> exponent_to_slot(
        static_cast<std::size_t>(*std::max_element(exponent_values.begin(), exponent_values.end()) + 1),
        -1);
    std::vector<std::vector<u64>> power_table;
    power_table.reserve(exponent_values.size());

    for (std::size_t slot = 0; slot < exponent_values.size(); ++slot) {
        const int e = exponent_values[slot];
        exponent_to_slot[static_cast<std::size_t>(e)] = static_cast<int>(slot);
        std::vector<u64> values(split_primes.size(), 0ULL);
        for (std::size_t i = 0; i < split_primes.size(); ++i) {
            const u64 pe = pow_with_limit(static_cast<u64>(split_primes[i]), e, limit);
            if (pe > limit) {
                break;
            }
            values[i] = pe;
        }
        power_table.push_back(std::move(values));
    }

    Evaluator evaluator;
    evaluator.n_limit = limit;
    evaluator.split_primes = split_primes;
    evaluator.exponent_values = exponent_values;
    evaluator.power_table = std::move(power_table);
    evaluator.exponent_to_slot = std::move(exponent_to_slot);
    evaluator.multiplier_table = multiplier_table;

    std::vector<Task> tasks;
    tasks.reserve(64);

    for (const auto& pattern : setup.exponent_patterns) {
        std::vector<int> perm = pattern;
        std::sort(perm.begin(), perm.end());
        do {
            if (perm.size() > split_primes.size()) {
                continue;
            }

            u64 product = 1ULL;
            bool feasible = true;
            for (std::size_t i = 0; i < perm.size(); ++i) {
                const auto& powv = evaluator.powers_for_exponent(perm[i]);
                const u64 pe = powv[i];
                if (pe == 0ULL || pe > limit / product) {
                    feasible = false;
                    break;
                }
                product *= pe;
            }

            if (feasible) {
                tasks.push_back(Task{perm});
            }
        } while (std::next_permutation(perm.begin(), perm.end()));
    }

    if (tasks.empty()) {
        return 0ULL;
    }

    std::vector<std::size_t> exponent_head_work(
        static_cast<std::size_t>(*std::max_element(exponent_values.begin(), exponent_values.end()) + 1),
        0ULL);
    for (const int e : exponent_values) {
        const auto& powv = evaluator.powers_for_exponent(e);
        std::size_t cnt = 0;
        while (cnt < powv.size() && powv[cnt] != 0ULL) {
            ++cnt;
        }
        exponent_head_work[static_cast<std::size_t>(e)] = cnt;
    }

    std::size_t workload_hint = 0ULL;
    for (const auto& task : tasks) {
        workload_hint += exponent_head_work[static_cast<std::size_t>(task.permutation.front())];
    }

    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, workload_hint);
    if (threads == 1U) {
        u64 total = 0ULL;
        for (const auto& task : tasks) {
            total += evaluator.evaluate_task(task);
        }
        return total;
    }

    std::atomic<std::size_t> next_task(0ULL);
    std::vector<std::thread> workers;
    std::vector<u64> partial(threads, 0ULL);
    workers.reserve(threads);

    for (unsigned t = 0; t < threads; ++t) {
        workers.emplace_back([&, t]() {
            u64 local = 0ULL;
            while (true) {
                const std::size_t idx = next_task.fetch_add(1ULL, std::memory_order_relaxed);
                if (idx >= tasks.size()) {
                    break;
                }
                local += evaluator.evaluate_task(tasks[idx]);
            }
            partial[t] = local;
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    u64 total = 0ULL;
    for (const u64 v : partial) {
        total += v;
    }
    return total;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints) {
        const u64 c1 = solve_f(kCheckpointN1,
                               kCheckpointR1,
                               options.allow_multithreading,
                               options.requested_threads);
        if (c1 != kCheckpointExpected1) {
            std::cerr << "Checkpoint failed: f(10^5, 4) = " << c1
                      << ", expected " << kCheckpointExpected1 << '\n';
            return 1;
        }

        const u64 c2 = solve_f(kCheckpointN2,
                               kCheckpointR2,
                               options.allow_multithreading,
                               options.requested_threads);
        if (c2 != kCheckpointExpected2) {
            std::cerr << "Checkpoint failed: f(10^8, 6) = " << c2
                      << ", expected " << kCheckpointExpected2 << '\n';
            return 1;
        }
    }

    const u64 answer = solve_f(options.limit,
                               options.r,
                               options.allow_multithreading,
                               options.requested_threads);
    std::cout << answer << '\n';
    return 0;
}
