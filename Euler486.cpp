#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;

constexpr u64 kMainMod = 87'654'321ULL;
constexpr u64 kDefaultLimit = 1'000'000'000'000'000'000ULL;
constexpr u64 kCheckpointL1 = 10'000'000ULL;
constexpr u64 kCheckpointL2 = 5'000'000'000ULL;

// For n = 9 + 6t + u (u in [0,5]):
// F5(n) = 2^(10+u) * 64^t - (200t + C[u]).
constexpr std::array<u64, 6> kC = {182ULL, 216ULL, 252ULL, 286ULL, 318ULL, 350ULL};
constexpr std::array<u64, 3> kPrimePowerFactors = {9ULL, 1997ULL, 4877ULL};

struct Options {
    u64 limit = kDefaultLimit;
    bool allow_multithreading = true;
    bool run_checkpoints = true;
    unsigned requested_threads = 0;
};

struct ModData {
    u64 mod = 1;
    u64 period = 1;
    std::array<std::vector<u32>, 6> residues_by_u;
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

        u64 limit = 0;
        if (parse_u64_after_prefix(arg, "--limit=", limit)) {
            options.limit = limit;
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

    return true;
}

u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((static_cast<u128>(a) * static_cast<u128>(b)) % mod);
}

u64 pow_mod(u64 base, u64 exp, u64 mod) {
    u64 result = 1 % mod;
    u64 cur = base % mod;
    u64 e = exp;

    while (e > 0) {
        if (e & 1ULL) {
            result = mul_mod(result, cur, mod);
        }
        cur = mul_mod(cur, cur, mod);
        e >>= 1ULL;
    }

    return result;
}

u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        const u64 t = a % b;
        a = b;
        b = t;
    }
    return a;
}

u64 lcm_u64(u64 a, u64 b) {
    return (a / gcd_u64(a, b)) * b;
}

std::vector<u64> distinct_prime_factors(u64 n) {
    std::vector<u64> factors;
    for (u64 p = 2; p * p <= n; ++p) {
        if (n % p != 0) {
            continue;
        }
        factors.push_back(p);
        while (n % p == 0) {
            n /= p;
        }
    }
    if (n > 1) {
        factors.push_back(n);
    }
    return factors;
}

u64 euler_phi(u64 n) {
    if (n == 0) {
        return 0;
    }

    u64 result = n;
    u64 x = n;
    for (u64 p = 2; p * p <= x; ++p) {
        if (x % p != 0) {
            continue;
        }
        while (x % p == 0) {
            x /= p;
        }
        result -= result / p;
    }
    if (x > 1) {
        result -= result / x;
    }
    return result;
}

u64 multiplicative_order(u64 a, u64 mod) {
    if (gcd_u64(a, mod) != 1ULL) {
        return 0;
    }

    u64 order = euler_phi(mod);
    const std::vector<u64> factors = distinct_prime_factors(order);
    for (const u64 p : factors) {
        while (order % p == 0ULL && pow_mod(a, order / p, mod) == 1ULL) {
            order /= p;
        }
    }
    return order;
}

unsigned choose_thread_count(bool allow_multithreading,
                             unsigned requested_threads,
                             std::size_t workload) {
    if (!allow_multithreading || workload < 2) {
        return 1;
    }

    unsigned threads = requested_threads;
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) {
            threads = 1;
        }
    }

    threads = std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(workload)));
    return threads;
}

// We keep enough prefix values to validate all F5 sample checkpoints.
std::vector<u64> compute_F5_prefix(int max_n) {
    std::vector<u64> f(static_cast<std::size_t>(max_n) + 1ULL, 0ULL);

    // State is last up to 5 bits.
    std::vector<std::pair<u32, u64>> states;
    states.emplace_back(0U, 1ULL);

    f[0] = 0ULL;

    for (int n = 1; n <= max_n; ++n) {
        u64 next_count[32] = {0ULL};
        for (const auto& [bits, count] : states) {
            for (u32 add = 0; add <= 1; ++add) {
                const u32 combined = (bits << 1U) | add;
                bool ok = true;
                if (n >= 5) {
                    const u32 v5 = combined & 31U;
                    const u32 b0 = (v5 >> 0U) & 1U;
                    const u32 b1 = (v5 >> 1U) & 1U;
                    const u32 b3 = (v5 >> 3U) & 1U;
                    const u32 b4 = (v5 >> 4U) & 1U;
                    if (b0 == b4 && b1 == b3) {
                        ok = false;
                    }
                }

                if (ok && n >= 6) {
                    const u32 v6 = combined & 63U;
                    const u32 b0 = (v6 >> 0U) & 1U;
                    const u32 b1 = (v6 >> 1U) & 1U;
                    const u32 b2 = (v6 >> 2U) & 1U;
                    const u32 b3 = (v6 >> 3U) & 1U;
                    const u32 b4 = (v6 >> 4U) & 1U;
                    const u32 b5 = (v6 >> 5U) & 1U;
                    if (b0 == b5 && b1 == b4 && b2 == b3) {
                        ok = false;
                    }
                }

                if (!ok) {
                    continue;
                }

                const u32 next_bits = (n <= 5) ? (combined & ((1U << n) - 1U)) : (combined & 31U);
                next_count[next_bits] += count;
            }
        }

        states.clear();
        states.reserve(32);
        u64 total_a = 0ULL;

        if (n <= 5) {
            const u32 limit = 1U << n;
            for (u32 mask = 0; mask < limit; ++mask) {
                if (next_count[mask] != 0ULL) {
                    states.emplace_back(mask, next_count[mask]);
                    total_a += next_count[mask];
                }
            }
        } else {
            for (u32 mask = 0; mask < 32U; ++mask) {
                if (next_count[mask] != 0ULL) {
                    states.emplace_back(mask, next_count[mask]);
                    total_a += next_count[mask];
                }
            }
        }

        const u64 total_strings = (n < 63) ? (1ULL << n) : 0ULL;
        f[static_cast<std::size_t>(n)] = f[static_cast<std::size_t>(n - 1)] + (total_strings - total_a);
    }

    return f;
}

std::vector<u32> compute_residues_for_mod_u(u64 mod, u64 period, int u) {
    std::vector<u32> residues;

    const u64 lhs_mul = pow_mod(2ULL, static_cast<u64>(10 + u), mod);
    u64 pow64 = 1ULL % mod;
    u64 rhs = kC[static_cast<std::size_t>(u)] % mod;

    residues.reserve(static_cast<std::size_t>(period / mod + 8ULL));

    for (u64 t = 0; t < period; ++t) {
        if (mul_mod(lhs_mul, pow64, mod) == rhs) {
            residues.push_back(static_cast<u32>(t));
        }

        pow64 = mul_mod(pow64, 64ULL % mod, mod);
        rhs += 200ULL;
        rhs %= mod;
    }

    return residues;
}

u64 mod_inverse(u64 a, u64 mod) {
    i64 t = 0;
    i64 new_t = 1;
    i64 r = static_cast<i64>(mod);
    i64 new_r = static_cast<i64>(a % mod);

    while (new_r != 0) {
        const i64 q = r / new_r;

        const i64 temp_t = t - q * new_t;
        t = new_t;
        new_t = temp_t;

        const i64 temp_r = r - q * new_r;
        r = new_r;
        new_r = temp_r;
    }

    if (r != 1) {
        return 0;
    }

    i64 x = t;
    x %= static_cast<i64>(mod);
    if (x < 0) {
        x += static_cast<i64>(mod);
    }
    return static_cast<u64>(x);
}

struct CRTPrecompute {
    u64 m1 = 1;
    u64 m2 = 1;
    u64 g = 1;
    u64 m1_div_g = 1;
    u64 m2_div_g = 1;
    u64 lcm = 1;
    u64 inv_m1_div_g_mod_m2_div_g = 0;
};

CRTPrecompute build_crt_precompute(u64 m1, u64 m2) {
    CRTPrecompute data;
    data.m1 = m1;
    data.m2 = m2;
    data.g = gcd_u64(m1, m2);
    data.m1_div_g = m1 / data.g;
    data.m2_div_g = m2 / data.g;
    data.lcm = data.m1_div_g * m2;
    data.inv_m1_div_g_mod_m2_div_g = mod_inverse(data.m1_div_g % data.m2_div_g, data.m2_div_g);
    return data;
}

u64 crt_merge_pair(u64 a, u64 b, const CRTPrecompute& crt) {
    // Assumes compatibility: (a - b) % g == 0.
    const i64 diff = static_cast<i64>(b) - static_cast<i64>(a);
    const i64 reduced = diff / static_cast<i64>(crt.g);

    i64 k = reduced % static_cast<i64>(crt.m2_div_g);
    if (k < 0) {
        k += static_cast<i64>(crt.m2_div_g);
    }

    const u64 ku = static_cast<u64>(k);
    const u64 t = mul_mod(ku, crt.inv_m1_div_g_mod_m2_div_g, crt.m2_div_g);
    const u64 x = static_cast<u64>((static_cast<u128>(a) +
                                    static_cast<u128>(crt.m1) * static_cast<u128>(t)) %
                                   crt.lcm);
    return x;
}

std::vector<u64> combine_residues_for_u(const std::vector<u32>& mod1997,
                                        u64 period1997,
                                        const std::vector<u32>& mod4877,
                                        u64 period4877,
                                        const std::vector<u32>& mod9,
                                        u64 period9,
                                        u64& out_global_period) {
    const CRTPrecompute crt12 = build_crt_precompute(period1997, period4877);
    const CRTPrecompute crt129 = build_crt_precompute(crt12.lcm, period9);

    out_global_period = crt129.lcm;

    std::vector<u32> mod4877_even;
    std::vector<u32> mod4877_odd;
    mod4877_even.reserve(mod4877.size() / 2 + 1);
    mod4877_odd.reserve(mod4877.size() / 2 + 1);
    for (const u32 value : mod4877) {
        if ((value & 1U) == 0U) {
            mod4877_even.push_back(value);
        } else {
            mod4877_odd.push_back(value);
        }
    }

    std::vector<u64> combined;
    combined.reserve(mod1997.size() * (mod4877.size() / 2 + 1) * std::max<std::size_t>(1, mod9.size()));

    for (const u32 a32 : mod1997) {
        const u64 a = static_cast<u64>(a32);
        const std::vector<u32>& candidates = ((a32 & 1U) == 0U) ? mod4877_even : mod4877_odd;

        for (const u32 b32 : candidates) {
            const u64 b = static_cast<u64>(b32);
            const u64 merged12 = crt_merge_pair(a, b, crt12);

            for (const u32 c32 : mod9) {
                const u64 c = static_cast<u64>(c32);
                const u64 merged = crt_merge_pair(merged12, c, crt129);
                combined.push_back(merged);
            }
        }
    }

    std::sort(combined.begin(), combined.end());
    combined.erase(std::unique(combined.begin(), combined.end()), combined.end());
    return combined;
}

u64 count_from_residues(const std::vector<u64>& residues, u64 period, u64 t_limit) {
    if (residues.empty()) {
        return 0ULL;
    }

    const u64 full = t_limit / period;
    const u64 rem = t_limit % period;

    const auto it = std::upper_bound(residues.begin(), residues.end(), rem);
    const u64 in_tail = static_cast<u64>(std::distance(residues.begin(), it));

    return full * static_cast<u64>(residues.size()) + in_tail;
}

class Euler486Solver {
public:
    explicit Euler486Solver(const Options& options)
        : allow_multithreading_(options.allow_multithreading),
          requested_threads_(options.requested_threads) {
        initialize_mod_data();
        build_mod_residue_tables();
        build_global_residue_tables();
    }

    u64 D(u64 limit) const {
        if (limit < 5ULL) {
            return 0ULL;
        }

        u64 answer = 0ULL;

        const u64 small_end = std::min<u64>(limit, 8ULL);
        for (u64 n = 5ULL; n <= small_end; ++n) {
            if (small_f5_mod_[static_cast<std::size_t>(n)] == 0ULL) {
                ++answer;
            }
        }

        if (limit < 9ULL) {
            return answer;
        }

        for (int u = 0; u < 6; ++u) {
            const u64 n0 = 9ULL + static_cast<u64>(u);
            if (n0 > limit) {
                continue;
            }

            const u64 t_limit = (limit - n0) / 6ULL;
            answer += count_from_residues(global_residues_by_u_[static_cast<std::size_t>(u)],
                                          global_period_,
                                          t_limit);
        }

        return answer;
    }

    bool run_checkpoints() const {
        bool ok = true;

        // F5 sample values.
        ok &= checkpoint_equal("F5(4)", small_f5_[4], 0ULL);
        ok &= checkpoint_equal("F5(5)", small_f5_[5], 8ULL);
        ok &= checkpoint_equal("F5(6)", small_f5_[6], 42ULL);
        ok &= checkpoint_equal("F5(11)", small_f5_[11], 3'844ULL);

        // D(L) sample values.
        ok &= checkpoint_equal("D(10^7)", D(kCheckpointL1), 0ULL);
        ok &= checkpoint_equal("D(5*10^9)", D(kCheckpointL2), 51ULL);

        // Additional structural checkpoint used by the derivation.
        static constexpr std::array<u64, 6> expected_a_cycle = {32ULL, 34ULL, 36ULL, 34ULL, 32ULL, 32ULL};
        for (int i = 0; i < 6; ++i) {
            const u64 n = 9ULL + static_cast<u64>(i);
            ok &= checkpoint_equal((std::string("A(") + std::to_string(n) + ")").c_str(),
                                  small_a_[static_cast<std::size_t>(n)],
                                  expected_a_cycle[static_cast<std::size_t>(i)]);
        }

        return ok;
    }

private:
    bool allow_multithreading_ = true;
    unsigned requested_threads_ = 0;

    std::array<ModData, 3> mod_data_;
    std::array<std::vector<u64>, 6> global_residues_by_u_;
    u64 global_period_ = 1ULL;

    std::vector<u64> small_f5_;
    std::vector<u64> small_a_;
    std::vector<u64> small_f5_mod_;

    static bool checkpoint_equal(const char* label, u64 got, u64 expected) {
        const bool ok = (got == expected);
        std::cout << "[checkpoint] " << label << " = " << got;
        if (!ok) {
            std::cout << " (expected " << expected << ")";
        }
        std::cout << '\n';
        return ok;
    }

    void initialize_mod_data() {
        for (std::size_t i = 0; i < mod_data_.size(); ++i) {
            mod_data_[i].mod = kPrimePowerFactors[i];
            const u64 order = multiplicative_order(64ULL % mod_data_[i].mod, mod_data_[i].mod);
            mod_data_[i].period = lcm_u64(mod_data_[i].mod, order);
        }

        small_f5_ = compute_F5_prefix(20);
        small_f5_mod_.assign(small_f5_.size(), 0ULL);
        for (std::size_t i = 0; i < small_f5_.size(); ++i) {
            small_f5_mod_[i] = small_f5_[i] % kMainMod;
        }

        // Also keep A(n) values for validation checks.
        small_a_.assign(21, 0ULL);
        small_a_[0] = 1ULL;
        for (std::size_t n = 1; n <= 20; ++n) {
            const u64 total = (n < 63) ? (1ULL << n) : 0ULL;
            small_a_[n] = total - (small_f5_[n] - small_f5_[n - 1]);
        }
    }

    void build_mod_residue_tables() {
        struct Task {
            int mod_index = 0;
            int u = 0;
        };

        std::vector<Task> tasks;
        tasks.reserve(18);
        for (int mi = 0; mi < 3; ++mi) {
            for (int u = 0; u < 6; ++u) {
                tasks.push_back(Task{mi, u});
            }
        }

        const unsigned threads = choose_thread_count(allow_multithreading_,
                                                     requested_threads_,
                                                     tasks.size());

        if (threads == 1U) {
            for (const Task& task : tasks) {
                ModData& md = mod_data_[static_cast<std::size_t>(task.mod_index)];
                md.residues_by_u[static_cast<std::size_t>(task.u)] =
                    compute_residues_for_mod_u(md.mod, md.period, task.u);
            }
            return;
        }

        std::atomic<std::size_t> next_task(0);
        std::vector<std::thread> workers;
        workers.reserve(threads);

        for (unsigned t = 0; t < threads; ++t) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next_task.fetch_add(1);
                    if (index >= tasks.size()) {
                        break;
                    }

                    const Task task = tasks[index];
                    ModData& md = mod_data_[static_cast<std::size_t>(task.mod_index)];
                    md.residues_by_u[static_cast<std::size_t>(task.u)] =
                        compute_residues_for_mod_u(md.mod, md.period, task.u);
                }
            });
        }

        for (auto& worker : workers) {
            worker.join();
        }
    }

    void build_global_residue_tables() {
        // Index mapping by modulus value.
        const int idx9 = 0;
        const int idx1997 = 1;
        const int idx4877 = 2;

        std::array<std::pair<int, std::vector<u64>>, 6> temp;
        for (int u = 0; u < 6; ++u) {
            temp[static_cast<std::size_t>(u)].first = u;
        }

        const unsigned threads = choose_thread_count(allow_multithreading_, requested_threads_, 6);
        if (threads == 1U) {
            for (int u = 0; u < 6; ++u) {
                u64 period = 1;
                temp[static_cast<std::size_t>(u)].second = combine_residues_for_u(
                    mod_data_[idx1997].residues_by_u[static_cast<std::size_t>(u)],
                    mod_data_[idx1997].period,
                    mod_data_[idx4877].residues_by_u[static_cast<std::size_t>(u)],
                    mod_data_[idx4877].period,
                    mod_data_[idx9].residues_by_u[static_cast<std::size_t>(u)],
                    mod_data_[idx9].period,
                    period);
                global_period_ = period;
            }
        } else {
            std::atomic<int> next_u(0);
            std::atomic<u64> period_seen(0ULL);
            std::vector<std::thread> workers;
            workers.reserve(threads);

            for (unsigned t = 0; t < threads; ++t) {
                workers.emplace_back([&]() {
                    while (true) {
                        const int u = next_u.fetch_add(1);
                        if (u >= 6) {
                            break;
                        }

                        u64 period = 1ULL;
                        auto combined = combine_residues_for_u(
                            mod_data_[idx1997].residues_by_u[static_cast<std::size_t>(u)],
                            mod_data_[idx1997].period,
                            mod_data_[idx4877].residues_by_u[static_cast<std::size_t>(u)],
                            mod_data_[idx4877].period,
                            mod_data_[idx9].residues_by_u[static_cast<std::size_t>(u)],
                            mod_data_[idx9].period,
                            period);

                        temp[static_cast<std::size_t>(u)].second = std::move(combined);

                        u64 expected = 0ULL;
                        if (period_seen.compare_exchange_strong(expected, period)) {
                            // first writer sets the period
                        } else {
                            if (expected != period) {
                                std::cerr << "Inconsistent global period across residues.\n";
                                std::terminate();
                            }
                        }
                    }
                });
            }

            for (auto& worker : workers) {
                worker.join();
            }

            global_period_ = period_seen.load();
        }

        for (const auto& item : temp) {
            global_residues_by_u_[static_cast<std::size_t>(item.first)] = item.second;
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Options options;
    if (!parse_arguments(argc, argv, options)) {
        std::cerr << "Usage: ./Euler486 [--limit=N] [--threads=T] [--single-thread] "
                     "[--skip-checkpoints]\n";
        return 1;
    }

    Euler486Solver solver(options);

    bool ok = true;
    if (options.run_checkpoints) {
        ok = solver.run_checkpoints();
    }

    const u64 answer = solver.D(options.limit);
    std::cout << "D(" << options.limit << ") = " << answer << '\n';

    if (!ok) {
        std::cerr << "Checkpoint validation failed.\n";
        return 2;
    }

    return 0;
}
