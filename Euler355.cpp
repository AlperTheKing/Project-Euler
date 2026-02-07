#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {

using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr u64 kDefaultN = 200'000ULL;

struct Options {
    u64 n = kDefaultN;
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
        if (parse_u64_after_prefix(arg, "--n=", parsed_u64)) {
            options.n = parsed_u64;
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
    if (!allow_multithreading || workload_units < 2ULL) {
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

std::vector<int> sieve_primes(const int n) {
    if (n < 2) {
        return {};
    }

    std::vector<bool> is_prime(static_cast<std::size_t>(n) + 1ULL, true);
    is_prime[0] = false;
    is_prime[1] = false;

    for (int i = 2; static_cast<u64>(i) * static_cast<u64>(i) <= static_cast<u64>(n); ++i) {
        if (!is_prime[static_cast<std::size_t>(i)]) {
            continue;
        }
        for (int j = i * i; j <= n; j += i) {
            is_prime[static_cast<std::size_t>(j)] = false;
        }
    }

    std::vector<int> primes;
    for (int i = 2; i <= n; ++i) {
        if (is_prime[static_cast<std::size_t>(i)]) {
            primes.push_back(i);
        }
    }
    return primes;
}

u64 max_prime_power_leq_n(const u64 n, const int prime) {
    u64 value = static_cast<u64>(prime);
    while (value <= n / static_cast<u64>(prime)) {
        value *= static_cast<u64>(prime);
    }
    return value;
}

u64 highest_power_with_limit(const int prime, const u64 limit) {
    if (limit < static_cast<u64>(prime)) {
        return 0ULL;
    }

    u64 value = static_cast<u64>(prime);
    while (value <= limit / static_cast<u64>(prime)) {
        value *= static_cast<u64>(prime);
    }
    return value;
}

struct SolverData {
    std::vector<int> all_primes;
    std::vector<u64> max_prime_powers;
    std::vector<int> small_primes;
    std::vector<u64> small_prime_powers;
    std::vector<int> large_primes;
    std::vector<u64> large_prime_powers;
    u64 baseline_sum = 1ULL;
};

SolverData prepare_solver_data(const u64 n) {
    SolverData data;
    data.all_primes = sieve_primes(static_cast<int>(n));

    data.max_prime_powers.reserve(data.all_primes.size());
    for (const int p : data.all_primes) {
        const u64 power = max_prime_power_leq_n(n, p);
        data.max_prime_powers.push_back(power);
        data.baseline_sum += power;
    }

    const int split = static_cast<int>(std::sqrt(static_cast<long double>(n)));

    for (std::size_t i = 0; i < data.all_primes.size(); ++i) {
        const int p = data.all_primes[i];
        const u64 power = data.max_prime_powers[i];
        if (p <= split) {
            data.small_primes.push_back(p);
            data.small_prime_powers.push_back(power);
        } else {
            data.large_primes.push_back(p);
            data.large_prime_powers.push_back(power);
        }
    }

    return data;
}

std::vector<std::tuple<int, int, u64>> build_gain_edges_chunk(const SolverData& data,
                                                               const u64 n,
                                                               const std::size_t begin,
                                                               const std::size_t end) {
    std::vector<std::tuple<int, int, u64>> edges;

    for (std::size_t si = begin; si < end; ++si) {
        const int p = data.small_primes[si];
        const u64 p_power = data.small_prime_powers[si];

        for (std::size_t li = 0; li < data.large_primes.size(); ++li) {
            const int q = data.large_primes[li];
            const u64 q_power = data.large_prime_powers[li];

            if (static_cast<u64>(p) * static_cast<u64>(q) > n) {
                break;
            }

            const u64 best_p_factor = highest_power_with_limit(p, n / static_cast<u64>(q));
            if (best_p_factor == 0ULL) {
                continue;
            }

            const u64 pair_value = best_p_factor * static_cast<u64>(q);
            if (pair_value <= p_power + q_power) {
                continue;
            }

            const u64 gain = pair_value - p_power - q_power;
            edges.emplace_back(static_cast<int>(si), static_cast<int>(li), gain);
        }
    }

    return edges;
}

std::vector<std::tuple<int, int, u64>> build_gain_edges(const SolverData& data,
                                                        const u64 n,
                                                        const bool allow_multithreading,
                                                        const unsigned requested_threads) {
    const std::size_t small_count = data.small_primes.size();
    if (small_count == 0ULL || data.large_primes.empty()) {
        return {};
    }

    const unsigned threads = choose_thread_count(allow_multithreading, requested_threads, small_count);
    if (threads == 1U) {
        return build_gain_edges_chunk(data, n, 0ULL, small_count);
    }

    std::vector<std::vector<std::tuple<int, int, u64>>> partial(
        static_cast<std::size_t>(threads));
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(threads));

    const std::size_t chunk = (small_count + static_cast<std::size_t>(threads) - 1ULL) /
                              static_cast<std::size_t>(threads);

    for (unsigned t = 0U; t < threads; ++t) {
        const std::size_t begin = static_cast<std::size_t>(t) * chunk;
        const std::size_t end = std::min<std::size_t>(small_count, begin + chunk);

        if (begin >= end) {
            break;
        }

        pool.emplace_back([&, t, begin, end]() {
            partial[static_cast<std::size_t>(t)] = build_gain_edges_chunk(data, n, begin, end);
        });
    }

    for (std::thread& th : pool) {
        th.join();
    }

    std::size_t total = 0ULL;
    for (const auto& part : partial) {
        total += part.size();
    }

    std::vector<std::tuple<int, int, u64>> edges;
    edges.reserve(total);
    for (auto& part : partial) {
        edges.insert(edges.end(), part.begin(), part.end());
    }

    return edges;
}

i64 hungarian_max_weight_bipartite(const int left_size,
                                   const int right_size,
                                   const std::vector<i64>& weights) {
    if (left_size == 0 || right_size == 0) {
        return 0;
    }

    const int n = left_size;
    const int m = right_size;

    auto cost = [&](const int i, const int j) -> i64 {
        const std::size_t idx = static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(m) +
                                static_cast<std::size_t>(j - 1);
        return -weights[idx];
    };

    const i64 inf = std::numeric_limits<i64>::max() / 4;

    std::vector<i64> u(static_cast<std::size_t>(n) + 1ULL, 0);
    std::vector<i64> v(static_cast<std::size_t>(m) + 1ULL, 0);
    std::vector<int> p(static_cast<std::size_t>(m) + 1ULL, 0);
    std::vector<int> way(static_cast<std::size_t>(m) + 1ULL, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;

        std::vector<i64> minv(static_cast<std::size_t>(m) + 1ULL, inf);
        std::vector<bool> used(static_cast<std::size_t>(m) + 1ULL, false);

        do {
            used[static_cast<std::size_t>(j0)] = true;
            const int i0 = p[static_cast<std::size_t>(j0)];

            i64 delta = inf;
            int j1 = 0;

            for (int j = 1; j <= m; ++j) {
                if (used[static_cast<std::size_t>(j)]) {
                    continue;
                }

                const i64 cur = cost(i0, j) - u[static_cast<std::size_t>(i0)] - v[static_cast<std::size_t>(j)];
                if (cur < minv[static_cast<std::size_t>(j)]) {
                    minv[static_cast<std::size_t>(j)] = cur;
                    way[static_cast<std::size_t>(j)] = j0;
                }

                if (minv[static_cast<std::size_t>(j)] < delta) {
                    delta = minv[static_cast<std::size_t>(j)];
                    j1 = j;
                }
            }

            for (int j = 0; j <= m; ++j) {
                if (used[static_cast<std::size_t>(j)]) {
                    u[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] += delta;
                    v[static_cast<std::size_t>(j)] -= delta;
                } else {
                    minv[static_cast<std::size_t>(j)] -= delta;
                }
            }

            j0 = j1;
        } while (p[static_cast<std::size_t>(j0)] != 0);

        do {
            const int j1 = way[static_cast<std::size_t>(j0)];
            p[static_cast<std::size_t>(j0)] = p[static_cast<std::size_t>(j1)];
            j0 = j1;
        } while (j0 != 0);
    }

    i64 min_cost = 0;
    for (int j = 1; j <= m; ++j) {
        if (p[static_cast<std::size_t>(j)] != 0) {
            min_cost += cost(p[static_cast<std::size_t>(j)], j);
        }
    }

    return -min_cost;
}

u64 solve_co_optimized(const u64 n,
                       const bool allow_multithreading,
                       const unsigned requested_threads) {
    const SolverData data = prepare_solver_data(n);
    const auto edges = build_gain_edges(data, n, allow_multithreading, requested_threads);

    const int left_size = static_cast<int>(data.small_primes.size());
    const int large_size = static_cast<int>(data.large_primes.size());
    const int right_size = large_size + left_size;  // includes dummy columns

    std::vector<i64> weights(static_cast<std::size_t>(left_size) * static_cast<std::size_t>(right_size), 0);

    for (const auto& edge : edges) {
        const int small_index = std::get<0>(edge);
        const int large_index = std::get<1>(edge);
        const i64 gain = static_cast<i64>(std::get<2>(edge));

        const std::size_t idx = static_cast<std::size_t>(small_index) * static_cast<std::size_t>(right_size) +
                                static_cast<std::size_t>(large_index);
        if (gain > weights[idx]) {
            weights[idx] = gain;
        }
    }

    const i64 max_gain = hungarian_max_weight_bipartite(left_size, right_size, weights);
    return data.baseline_sum + static_cast<u64>(max_gain);
}

bool run_checkpoints() {
    struct Checkpoint {
        u64 n;
        u64 expected;
    };

    const std::vector<Checkpoint> checkpoints = {
        {10ULL, 30ULL},
        {30ULL, 193ULL},
        {100ULL, 1356ULL},
    };

    for (const Checkpoint& cp : checkpoints) {
        const u64 got = solve_co_optimized(cp.n, false, 1U);
        if (got != cp.expected) {
            std::cerr << "Checkpoint failed for n=" << cp.n << ": expected " << cp.expected
                      << ", got " << got << '\n';
            return false;
        }
    }

    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_arguments(argc, argv, options)) {
        return 1;
    }

    if (options.run_checkpoints && !run_checkpoints()) {
        return 1;
    }

    const u64 answer = solve_co_optimized(options.n,
                                          options.allow_multithreading,
                                          options.requested_threads);

    std::cout << answer << '\n';
    return 0;
}
